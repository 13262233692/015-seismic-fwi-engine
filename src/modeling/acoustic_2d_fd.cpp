#include "fwi/modeling/acoustic_2d_fd.h"
#include "fwi/common/utils.h"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace fwi {

static const float64 FD_COEFF_4_ORDER[5] = {
    0.0,
    9.0 / 8.0,
    -1.0 / 24.0,
    0.0,
    0.0
};

static const float64 FD_COEFF_4_ORDER_STAGGERED[4] = {
    1.0 / 2.0 + 1.0 / 24.0,
    -1.0 / 24.0,
    0.0,
    0.0
};

Acoustic2DFD::Acoustic2DFD()
    : model_()
    , params_()
    , nx_(0), nz_(0), nx_total_(0), nz_total_(0)
    , dx_(0.0), dz_(0.0)
    , pressure_(), pressure_prev_()
    , velocity_x_(), velocity_z_()
    , density_(), modulus_()
    , pml_damping_x_(), pml_damping_z_()
    , pml_psi_x_(), pml_psi_z_()
    , sources_(), receivers_(), receiver_indices_()
    , current_step_(0), current_time_(0.0)
    , snapshots_(), receiver_data_(), time_axis_()
    , snapshot_callback_(nullptr)
    , step_callback_(nullptr)
    , fd_order_(4)
{}

Acoustic2DFD::Acoustic2DFD(const VelocityModel2D& model, const FDSimulationParams& params)
    : Acoustic2DFD()
{
    model_ = model;
    params_ = params;
}

void Acoustic2DFD::set_model(const VelocityModel2D& model) {
    model_ = model;
}

void Acoustic2DFD::set_params(const FDSimulationParams& params) {
    params_ = params;
}

const VelocityModel2D& Acoustic2DFD::model() const { return model_; }
const FDSimulationParams& Acoustic2DFD::params() const { return params_; }

const std::vector<Snapshot2D>& Acoustic2DFD::snapshots() const { return snapshots_; }
const std::vector<std::vector<float64>>& Acoustic2DFD::receiver_data() const { return receiver_data_; }
const std::vector<float64>& Acoustic2DFD::time_axis() const { return time_axis_; }

void Acoustic2DFD::add_source(const Source& source) {
    sources_.push_back(source);
}

void Acoustic2DFD::set_receivers(const std::vector<Receiver>& receivers) {
    receivers_ = receivers;
}

const Array2D<float64>& Acoustic2DFD::current_pressure() const { return pressure_; }
const Array2D<float64>& Acoustic2DFD::current_velocity_x() const { return velocity_x_; }
const Array2D<float64>& Acoustic2DFD::current_velocity_z() const { return velocity_z_; }

float64 Acoustic2DFD::current_time() const { return current_time_; }
int32 Acoustic2DFD::current_time_step() const { return current_step_; }

void Acoustic2DFD::set_snapshot_callback(SnapshotCallback callback) {
    snapshot_callback_ = callback;
}

void Acoustic2DFD::set_step_callback(StepCallback callback) {
    step_callback_ = callback;
}

void Acoustic2DFD::init_pml_coefficients() {
    if (!params_.use_pml || params_.pml_width <= 0) {
        return;
    }

    int32 npml = params_.pml_width;
    float64 vmax = model_.get_max_velocity();
    float64 freq_max = 0.0;
    for (const auto& src : sources_) {
        freq_max = std::max(freq_max, src.wavelet().peak_frequency());
    }
    if (freq_max <= 0.0) freq_max = 50.0;

    float64 lambda_min = vmax / freq_max;
    float64 d0 = 3.0 * vmax / (2.0 * static_cast<float64>(npml) * std::min(dx_, dz_));

    pml_damping_x_.fill(0.0);
    pml_damping_z_.fill(0.0);

    for (int32 ix = 0; ix < nx_total_; ix++) {
        for (int32 iz = 0; iz < nz_total_; iz++) {
            float64 damp_x = 0.0;
            float64 damp_z = 0.0;

            if (ix < npml) {
                float64 norm_dist = static_cast<float64>(npml - ix) / npml;
                damp_x = d0 * norm_dist * norm_dist;
            } else if (ix >= nx_total_ - npml) {
                float64 norm_dist = static_cast<float64>(ix - (nx_total_ - npml - 1)) / npml;
                damp_x = d0 * norm_dist * norm_dist;
            }

            if (iz < npml) {
                float64 norm_dist = static_cast<float64>(npml - iz) / npml;
                damp_z = d0 * norm_dist * norm_dist;
            } else if (iz >= nz_total_ - npml) {
                float64 norm_dist = static_cast<float64>(iz - (nz_total_ - npml - 1)) / npml;
                damp_z = d0 * norm_dist * norm_dist;
            }

            pml_damping_x_(ix, iz) = damp_x;
            pml_damping_z_(ix, iz) = damp_z;
        }
    }
}

void Acoustic2DFD::apply_pml_damping(float64& val_x, float64& val_z, int32 ix, int32 iz) {
    if (!params_.use_pml || params_.pml_width <= 0) {
        return;
    }

    float64 dx = pml_damping_x_(ix, iz);
    float64 dz = pml_damping_z_(ix, iz);
    float64 dt = params_.dt;

    val_x *= std::exp(-dx * dt);
    val_z *= std::exp(-dz * dt);
}

void Acoustic2DFD::initialize() {
    const Grid2D& grid = model_.grid();
    nx_ = grid.nx;
    nz_ = grid.nz;
    dx_ = grid.dx;
    dz_ = grid.dz;

    if (params_.use_pml && params_.pml_width > 0) {
        nx_total_ = nx_ + 2 * params_.pml_width;
        nz_total_ = nz_ + 2 * params_.pml_width;
    } else {
        nx_total_ = nx_;
        nz_total_ = nz_;
    }

    pressure_ = Array2D<float64>(nx_total_, nz_total_, 0.0);
    pressure_prev_ = Array2D<float64>(nx_total_, nz_total_, 0.0);
    velocity_x_ = Array2D<float64>(nx_total_, nz_total_, 0.0);
    velocity_z_ = Array2D<float64>(nx_total_, nz_total_, 0.0);
    density_ = Array2D<float64>(nx_total_, nz_total_, 1000.0);
    modulus_ = Array2D<float64>(nx_total_, nz_total_, 0.0);

    if (params_.use_pml && params_.pml_width > 0) {
        pml_damping_x_ = Array2D<float64>(nx_total_, nz_total_, 0.0);
        pml_damping_z_ = Array2D<float64>(nx_total_, nz_total_, 0.0);
        pml_psi_x_ = Array2D<float64>(nx_total_, nz_total_, 0.0);
        pml_psi_z_ = Array2D<float64>(nx_total_, nz_total_, 0.0);
    }

    int32 pml = params_.use_pml ? params_.pml_width : 0;
    for (int32 ix = 0; ix < nx_; ix++) {
        for (int32 iz = 0; iz < nz_; iz++) {
            int32 ix_total = ix + pml;
            int32 iz_total = iz + pml;
            float64 vel = model_.get_velocity(ix, iz);
            float64 rho = model_.get_density(ix, iz);
            density_(ix_total, iz_total) = rho;
            modulus_(ix_total, iz_total) = rho * vel * vel;
        }
    }

    for (int32 ix = 0; ix < nx_total_; ix++) {
        for (int32 iz = 0; iz < nz_total_; iz++) {
            int32 ix_model = clamp(ix - pml, 0, nx_ - 1);
            int32 iz_model = clamp(iz - pml, 0, nz_ - 1);
            if (density_(ix, iz) <= 0.0) {
                density_(ix, iz) = model_.get_density(ix_model, iz_model);
            }
            if (modulus_(ix, iz) <= 0.0) {
                float64 vel = model_.get_velocity(ix_model, iz_model);
                modulus_(ix, iz) = density_(ix, iz) * vel * vel;
            }
        }
    }

    receiver_indices_.clear();
    Grid2D pml_grid(nx_total_, nz_total_, dx_, dz_,
                    grid.ox - pml * dx_, grid.oz - pml * dz_);
    for (const auto& recv : receivers_) {
        int32 ix = nearest_index(recv.x(), pml_grid.ox, pml_grid.dx);
        int32 iz = nearest_index(recv.z(), pml_grid.oz, pml_grid.dz);
        ix = clamp(ix, 0, nx_total_ - 1);
        iz = clamp(iz, 0, nz_total_ - 1);
        receiver_indices_.emplace_back(ix, iz);
    }

    if (params_.record_receivers) {
        receiver_data_.resize(receivers_.size());
        for (auto& data : receiver_data_) {
            data.reserve(params_.num_time_steps);
        }
    }

    time_axis_.clear();
    time_axis_.reserve(params_.num_time_steps);

    snapshots_.clear();

    if (params_.use_pml && params_.pml_width > 0) {
        init_pml_coefficients();
    }

    current_step_ = 0;
    current_time_ = 0.0;
}

void Acoustic2DFD::apply_source(int32 time_idx) {
    int32 pml = params_.use_pml ? params_.pml_width : 0;
    Grid2D pml_grid(nx_total_, nz_total_, dx_, dz_,
                    model_.grid().ox - pml * dx_, model_.grid().oz - pml * dz_);

    for (const auto& src : sources_) {
        float64 amplitude = src.at(time_idx);
        if (std::abs(amplitude) < 1e-15) continue;

        int32 ix, iz;
        src.get_grid_indices(pml_grid, ix, iz);

        switch (src.type()) {
            case SourceType::EXPLOSIVE:
            case SourceType::PRESSURE: {
                float64 source_term = amplitude * modulus_(ix, iz) * params_.dt / density_(ix, iz);
                pressure_(ix, iz) += source_term;
                break;
            }
            case SourceType::PARTICLE_VELOCITY_X: {
                velocity_x_(ix, iz) += amplitude * params_.dt;
                break;
            }
            case SourceType::PARTICLE_VELOCITY_Z: {
                velocity_z_(ix, iz) += amplitude * params_.dt;
                break;
            }
        }
    }
}

void Acoustic2DFD::record_receivers(int32 time_idx) {
    if (!params_.record_receivers || receiver_indices_.empty()) {
        return;
    }

    if (time_idx % params_.record_interval != 0) {
        return;
    }

    for (size_t i = 0; i < receiver_indices_.size(); i++) {
        int32 ix = receiver_indices_[i].first;
        int32 iz = receiver_indices_[i].second;
        receiver_data_[i].push_back(pressure_(ix, iz));
    }
}

void Acoustic2DFD::save_snapshot(int32 time_idx) {
    if (params_.snapshot_interval <= 0 || time_idx % params_.snapshot_interval != 0) {
        return;
    }

    int32 pml = params_.use_pml ? params_.pml_width : 0;

    Snapshot2D snap;
    snap.time_step = time_idx;
    snap.time = current_time_;
    snap.pressure = Array2D<float64>(nx_, nz_, 0.0);
    snap.velocity_x = Array2D<float64>(nx_, nz_, 0.0);
    snap.velocity_z = Array2D<float64>(nx_, nz_, 0.0);

    for (int32 ix = 0; ix < nx_; ix++) {
        for (int32 iz = 0; iz < nz_; iz++) {
            int32 ix_total = ix + pml;
            int32 iz_total = iz + pml;
            snap.pressure(ix, iz) = pressure_(ix_total, iz_total);
            snap.velocity_x(ix, iz) = velocity_x_(ix_total, iz_total);
            snap.velocity_z(ix, iz) = velocity_z_(ix_total, iz_total);
        }
    }

    snapshots_.push_back(snap);

    if (snapshot_callback_) {
        snapshot_callback_(snap);
    }
}

void Acoustic2DFD::update_velocity_x(int32 nx, int32 nz, Array2D<float64>& vx, const Array2D<float64>& p,
                                     const Array2D<float64>& density, float64 dt, float64 dx) {
    float64 coeff1 = FD_COEFF_4_ORDER_STAGGERED[0];
    float64 coeff2 = FD_COEFF_4_ORDER_STAGGERED[1];
    float64 dt_over_dx = dt / dx;

    for (int32 iz = 2; iz < nz - 2; iz++) {
        for (int32 ix = 2; ix < nx - 2; ix++) {
            float64 dp_dx = coeff1 * (p(ix, iz) - p(ix - 1, iz))
                          + coeff2 * (p(ix + 1, iz) - p(ix - 2, iz));
            dp_dx *= dt_over_dx;
            float64 rho = density(ix, iz);
            vx(ix, iz) -= dp_dx / rho;
        }
    }
}

void Acoustic2DFD::update_velocity_z(int32 nx, int32 nz, Array2D<float64>& vz, const Array2D<float64>& p,
                                     const Array2D<float64>& density, float64 dt, float64 dz) {
    float64 coeff1 = FD_COEFF_4_ORDER_STAGGERED[0];
    float64 coeff2 = FD_COEFF_4_ORDER_STAGGERED[1];
    float64 dt_over_dz = dt / dz;

    for (int32 iz = 2; iz < nz - 2; iz++) {
        for (int32 ix = 2; ix < nx - 2; ix++) {
            float64 dp_dz = coeff1 * (p(ix, iz) - p(ix, iz - 1))
                          + coeff2 * (p(ix, iz + 1) - p(ix, iz - 2));
            dp_dz *= dt_over_dz;
            float64 rho = density(ix, iz);
            vz(ix, iz) -= dp_dz / rho;
        }
    }
}

void Acoustic2DFD::update_pressure(int32 nx, int32 nz, Array2D<float64>& p,
                                   const Array2D<float64>& vx, const Array2D<float64>& vz,
                                   const Array2D<float64>& modulus, float64 dt, float64 dx, float64 dz) {
    float64 coeff1 = FD_COEFF_4_ORDER_STAGGERED[0];
    float64 coeff2 = FD_COEFF_4_ORDER_STAGGERED[1];
    float64 dt_over_dx = dt / dx;
    float64 dt_over_dz = dt / dz;

    for (int32 iz = 2; iz < nz - 2; iz++) {
        for (int32 ix = 2; ix < nx - 2; ix++) {
            float64 dvx_dx = coeff1 * (vx(ix + 1, iz) - vx(ix, iz))
                           + coeff2 * (vx(ix + 2, iz) - vx(ix - 1, iz));
            dvx_dx *= dt_over_dx;

            float64 dvz_dz = coeff1 * (vz(ix, iz + 1) - vz(ix, iz))
                           + coeff2 * (vz(ix, iz + 2) - vz(ix, iz - 1));
            dvz_dz *= dt_over_dz;

            float64 K = modulus(ix, iz);
            p(ix, iz) -= K * (dvx_dx + dvz_dz);
        }
    }
}

void Acoustic2DFD::step(int32 num_steps) {
    for (int32 step = 0; step < num_steps; step++) {
        if (current_step_ >= params_.num_time_steps) {
            return;
        }

        apply_source(current_step_);

        update_velocity_x(nx_total_, nz_total_, velocity_x_, pressure_, density_, params_.dt, dx_);
        update_velocity_z(nx_total_, nz_total_, velocity_z_, pressure_, density_, params_.dt, dz_);

        if (params_.use_pml && params_.pml_width > 0) {
            for (int32 iz = 0; iz < nz_total_; iz++) {
                for (int32 ix = 0; ix < nx_total_; ix++) {
                    float64 vx = velocity_x_(ix, iz);
                    float64 vz = velocity_z_(ix, iz);
                    apply_pml_damping(vx, vz, ix, iz);
                    velocity_x_(ix, iz) = vx;
                    velocity_z_(ix, iz) = vz;
                }
            }
        }

        update_pressure(nx_total_, nz_total_, pressure_, velocity_x_, velocity_z_, modulus_, params_.dt, dx_, dz_);

        if (params_.use_pml && params_.pml_width > 0) {
            for (int32 iz = 0; iz < nz_total_; iz++) {
                for (int32 ix = 0; ix < nx_total_; ix++) {
                    float64 p = pressure_(ix, iz);
                    float64 dummy = 0.0;
                    apply_pml_damping(p, dummy, ix, iz);
                    pressure_(ix, iz) = p;
                }
            }
        }

        current_step_++;
        current_time_ = current_step_ * params_.dt;
        time_axis_.push_back(current_time_);

        record_receivers(current_step_);
        save_snapshot(current_step_);

        if (step_callback_) {
            step_callback_(current_step_, current_time_);
        }
    }
}

void Acoustic2DFD::run() {
    initialize();

    save_snapshot(0);

    while (current_step_ < params_.num_time_steps) {
        step(1);
    }
}

void Acoustic2DFD::compute_derivatives_x(const Array2D<float64>& field, Array2D<float64>& deriv, int32 order) {
    int32 nx = field.nx();
    int32 nz = field.nz();
    float64 inv_2dx = 1.0 / (2.0 * dx_);

    for (int32 iz = 1; iz < nz - 1; iz++) {
        for (int32 ix = 1; ix < nx - 1; ix++) {
            deriv(ix, iz) = (field(ix + 1, iz) - field(ix - 1, iz)) * inv_2dx;
        }
    }
}

void Acoustic2DFD::compute_derivatives_z(const Array2D<float64>& field, Array2D<float64>& deriv, int32 order) {
    int32 nx = field.nx();
    int32 nz = field.nz();
    float64 inv_2dz = 1.0 / (2.0 * dz_);

    for (int32 iz = 1; iz < nz - 1; iz++) {
        for (int32 ix = 1; ix < nx - 1; ix++) {
            deriv(ix, iz) = (field(ix, iz + 1) - field(ix, iz - 1)) * inv_2dz;
        }
    }
}

}
