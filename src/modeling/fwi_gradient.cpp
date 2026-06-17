#include "fwi/modeling/fwi_gradient.h"
#include "fwi/common/utils.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <cstring>

namespace fwi {

FwiGradient::FwiGradient()
    : model_()
    , params_()
    , wavelet_()
    , geometry_()
    , num_threads_(4)
    , mode_(GradientComputationMode::OPENMP_REDUCTION)
    , last_gradient_()
    , last_residual_(0.0)
    , nx_(0), nz_(0), nx_total_(0), nz_total_(0), pml_(0)
    , dx_(0.0), dz_(0.0), dt_(0.0)
{
}

void FwiGradient::set_velocity_model(const VelocityModel2D& model) {
    model_ = model;
    nx_ = model.grid().nx;
    nz_ = model.grid().nz;
    dx_ = model.grid().dx;
    dz_ = model.grid().dz;
}

void FwiGradient::set_params(const FDSimulationParams& params) {
    params_ = params;
    dt_ = params.dt;
    pml_ = params.use_pml ? params.pml_width : 0;
    nx_total_ = nx_ + 2 * pml_;
    nz_total_ = nz_ + 2 * pml_;
}

void FwiGradient::set_wavelet(const Wavelet& wavelet) {
    wavelet_ = wavelet;
}

void FwiGradient::set_geometry(const AcquisitionGeometry& geometry) {
    geometry_ = geometry;
}

void FwiGradient::set_num_threads(int32 num_threads) {
    num_threads_ = num_threads > 0 ? num_threads : 1;
    omp_set_num_threads(num_threads_);
}

int32 FwiGradient::num_threads() const {
    return num_threads_;
}

void FwiGradient::set_computation_mode(GradientComputationMode mode) {
    mode_ = mode;
}

GradientComputationMode FwiGradient::computation_mode() const {
    return mode_;
}

const Array2D<float64>& FwiGradient::last_gradient() const {
    return last_gradient_;
}

double FwiGradient::last_residual() const {
    return last_residual_;
}

void FwiGradient::compute_residual(
    const std::vector<std::vector<float64>>& observed,
    const std::vector<std::vector<float64>>& synthetic,
    std::vector<std::vector<float64>>& residual,
    double& residual_norm
) {
    size_t nrecv = observed.size();
    residual.resize(nrecv);
    residual_norm = 0.0;

    for (size_t ir = 0; ir < nrecv; ir++) {
        size_t nt = observed[ir].size();
        residual[ir].resize(nt, 0.0);
        for (size_t it = 0; it < nt; it++) {
            double syn = (it < synthetic[ir].size()) ? synthetic[ir][it] : 0.0;
            residual[ir][it] = observed[ir][it] - syn;
            residual_norm += residual[ir][it] * residual[ir][it];
        }
    }
    residual_norm = std::sqrt(residual_norm);
}

std::vector<Snapshot2D> FwiGradient::run_forward(int32 shot_index) {
    Acoustic2DFD solver(model_, params_);

    const Shot& shot = geometry_.shot(shot_index);
    Source source(shot, wavelet_, SourceType::EXPLOSIVE);
    solver.add_source(source);
    solver.set_receivers(geometry_.receivers());

    FDSimulationParams save_params = params_;
    save_params.snapshot_interval = 1;
    solver.set_params(save_params);

    solver.run();

    return solver.snapshots();
}

std::vector<std::vector<float64>> FwiGradient::run_adjoint(
    int32 shot_index,
    const std::vector<std::vector<float64>>& residual,
    std::vector<Snapshot2D>& adjoint_wavefield_out
) {
    int32 nt = params_.num_time_steps;
    int32 nrecv = static_cast<int32>(geometry_.num_receivers());

    Acoustic2DFD solver(model_, params_);
    const Shot& shot = geometry_.shot(shot_index);

    for (int32 ir = 0; ir < nrecv; ir++) {
        const Receiver& recv = geometry_.receiver(ir);

        std::vector<float64> adjoint_src_values(nt, 0.0);
        for (int32 it = 0; it < nt; it++) {
            int32 t_rev = nt - 1 - it;
            if (ir < static_cast<int32>(residual.size()) &&
                t_rev < static_cast<int32>(residual[ir].size())) {
                adjoint_src_values[it] = residual[ir][t_rev];
            }
        }

        Wavelet recv_wavelet = wavelet_;
        recv_wavelet.data() = adjoint_src_values;

        Shot recv_as_shot(recv.id(), recv.x(), recv.y(), recv.z());
        Source adjoint_source(recv_as_shot, recv_wavelet, SourceType::EXPLOSIVE);
        solver.add_source(adjoint_source);
    }

    Receiver shot_as_receiver(shot.id(), shot.x(), shot.y(), shot.z());
    solver.set_receivers({shot_as_receiver});

    FDSimulationParams save_params = params_;
    save_params.snapshot_interval = 1;
    solver.set_params(save_params);

    solver.run();

    adjoint_wavefield_out = solver.snapshots();

    return solver.receiver_data();
}

void FwiGradient::compute_gradient_serial(
    const std::vector<Snapshot2D>& forward_wavefield,
    const std::vector<Snapshot2D>& adjoint_wavefield,
    Array2D<float64>& gradient
) {
    int32 nt = static_cast<int32>(forward_wavefield.size());
    int32 nt_adj = static_cast<int32>(adjoint_wavefield.size());

    gradient.fill(0.0);

    for (int32 t = 0; t < nt; t++) {
        int32 t_adj = nt_adj - 1 - t;
        if (t_adj < 0) t_adj = 0;

        const Array2D<float64>& p_fwd = forward_wavefield[t].pressure;
        const Array2D<float64>& p_adj = adjoint_wavefield[t_adj].pressure;

        for (int32 iz = 0; iz < nz_; iz++) {
            for (int32 ix = 0; ix < nx_; ix++) {
                double vel = model_.get_velocity(ix, iz);
                double v3 = vel * vel * vel;
                double coef = -2.0 * dt_ / v3;
                gradient(ix, iz) += coef * p_fwd(ix, iz) * p_adj(ix, iz);
            }
        }
    }
}

void FwiGradient::compute_gradient_openmp_unsafe(
    const std::vector<Snapshot2D>& forward_wavefield,
    const std::vector<Snapshot2D>& adjoint_wavefield,
    Array2D<float64>& gradient
) {
    int32 nt = static_cast<int32>(forward_wavefield.size());
    int32 nt_adj = static_cast<int32>(adjoint_wavefield.size());

    gradient.fill(0.0);

    #pragma omp parallel for num_threads(num_threads_) schedule(dynamic, 1)
    for (int32 t = 0; t < nt; t++) {
        int32 t_adj = nt_adj - 1 - t;
        if (t_adj < 0) t_adj = 0;

        const Array2D<float64>& p_fwd = forward_wavefield[t].pressure;
        const Array2D<float64>& p_adj = adjoint_wavefield[t_adj].pressure;

        for (int32 iz = 0; iz < nz_; iz++) {
            for (int32 ix = 0; ix < nx_; ix++) {
                double vel = model_.get_velocity(ix, iz);
                double v3 = vel * vel * vel;
                double coef = -2.0 * dt_ / v3;
                double contrib = coef * p_fwd(ix, iz) * p_adj(ix, iz);

                gradient(ix, iz) += contrib;
            }
        }
    }
}

void FwiGradient::compute_gradient_openmp_atomic(
    const std::vector<Snapshot2D>& forward_wavefield,
    const std::vector<Snapshot2D>& adjoint_wavefield,
    Array2D<float64>& gradient
) {
    int32 nt = static_cast<int32>(forward_wavefield.size());
    int32 nt_adj = static_cast<int32>(adjoint_wavefield.size());

    gradient.fill(0.0);

    #pragma omp parallel for num_threads(num_threads_)
    for (int32 t = 0; t < nt; t++) {
        int32 t_adj = nt_adj - 1 - t;
        if (t_adj < 0) t_adj = 0;

        const Array2D<float64>& p_fwd = forward_wavefield[t].pressure;
        const Array2D<float64>& p_adj = adjoint_wavefield[t_adj].pressure;

        for (int32 iz = 0; iz < nz_; iz++) {
            for (int32 ix = 0; ix < nx_; ix++) {
                double vel = model_.get_velocity(ix, iz);
                double v3 = vel * vel * vel;
                double coef = -2.0 * dt_ / v3;
                double contrib = coef * p_fwd(ix, iz) * p_adj(ix, iz);

                #pragma omp atomic
                gradient(ix, iz) += contrib;
            }
        }
    }
}

void FwiGradient::compute_gradient_openmp_reduction(
    const std::vector<Snapshot2D>& forward_wavefield,
    const std::vector<Snapshot2D>& adjoint_wavefield,
    Array2D<float64>& gradient
) {
    int32 nt = static_cast<int32>(forward_wavefield.size());
    int32 nt_adj = static_cast<int32>(adjoint_wavefield.size());
    int64_t ngrid = static_cast<int64_t>(nx_) * nz_;

    gradient.fill(0.0);
    double* grad_ptr = gradient.data();

    int32 actual_threads = num_threads_;
    std::vector<std::vector<double>> thread_local_grads(
        actual_threads, std::vector<double>(ngrid, 0.0)
    );

    #pragma omp parallel num_threads(actual_threads)
    {
        int32 tid = omp_get_thread_num();
        std::vector<double>& local_grad = thread_local_grads[tid];

        #pragma omp for schedule(static)
        for (int32 t = 0; t < nt; t++) {
            int32 t_adj = nt_adj - 1 - t;
            if (t_adj < 0) t_adj = 0;

            const Array2D<float64>& p_fwd = forward_wavefield[t].pressure;
            const Array2D<float64>& p_adj = adjoint_wavefield[t_adj].pressure;

            for (int32 iz = 0; iz < nz_; iz++) {
                for (int32 ix = 0; ix < nx_; ix++) {
                    double vel = model_.get_velocity(ix, iz);
                    double v3 = vel * vel * vel;
                    double coef = -2.0 * dt_ / v3;
                    int64_t idx = static_cast<int64_t>(iz) * nx_ + ix;
                    local_grad[idx] += coef * p_fwd(ix, iz) * p_adj(ix, iz);
                }
            }
        }
    }

    for (int32 tid = 0; tid < actual_threads; tid++) {
        const std::vector<double>& local_grad = thread_local_grads[tid];
        for (int64_t idx = 0; idx < ngrid; idx++) {
            grad_ptr[idx] += local_grad[idx];
        }
    }
}

FwiGradientResult FwiGradient::compute_gradient(
    int32 shot_index,
    const std::vector<std::vector<float64>>& observed_data,
    const std::vector<Snapshot2D>& forward_wavefield_in
) {
    if (params_.dt <= 0.0) {
        params_.dt = model_.compute_stable_dt(0.5);
        dt_ = params_.dt;
    }
    if (pml_ > 0) {
        nx_total_ = nx_ + 2 * pml_;
        nz_total_ = nz_ + 2 * pml_;
    } else {
        nx_total_ = nx_;
        nz_total_ = nz_;
    }

    std::vector<Snapshot2D> forward_wavefield;
    std::vector<std::vector<float64>> synthetic_data;

    if (forward_wavefield_in.empty()) {
        Acoustic2DFD solver(model_, params_);
        const Shot& shot = geometry_.shot(shot_index);
        Source source(shot, wavelet_, SourceType::EXPLOSIVE);
        solver.add_source(source);
        solver.set_receivers(geometry_.receivers());

        FDSimulationParams save_params = params_;
        save_params.snapshot_interval = 1;
        solver.set_params(save_params);

        solver.run();

        forward_wavefield = solver.snapshots();
        synthetic_data = solver.receiver_data();
    } else {
        forward_wavefield = forward_wavefield_in;
        Acoustic2DFD solver(model_, params_);
        const Shot& shot = geometry_.shot(shot_index);
        Source source(shot, wavelet_, SourceType::EXPLOSIVE);
        solver.add_source(source);
        solver.set_receivers(geometry_.receivers());
        solver.run();
        synthetic_data = solver.receiver_data();
    }

    std::vector<std::vector<float64>> residual;
    compute_residual(observed_data, synthetic_data, residual, last_residual_);

    std::vector<Snapshot2D> adjoint_wavefield;
    run_adjoint(shot_index, residual, adjoint_wavefield);

    last_gradient_ = Array2D<float64>(nx_, nz_, 0.0);

    switch (mode_) {
        case GradientComputationMode::SERIAL:
            compute_gradient_serial(forward_wavefield, adjoint_wavefield, last_gradient_);
            break;
        case GradientComputationMode::OPENMP_UNSAFE:
            compute_gradient_openmp_unsafe(forward_wavefield, adjoint_wavefield, last_gradient_);
            break;
        case GradientComputationMode::OPENMP_ATOMIC:
            compute_gradient_openmp_atomic(forward_wavefield, adjoint_wavefield, last_gradient_);
            break;
        case GradientComputationMode::OPENMP_REDUCTION:
            compute_gradient_openmp_reduction(forward_wavefield, adjoint_wavefield, last_gradient_);
            break;
    }

    FwiGradientResult result;
    result.velocity_gradient = last_gradient_;
    result.residual_norm = last_residual_;
    result.num_time_steps = params_.num_time_steps;

    return result;
}

FwiGradientResult FwiGradient::compute_gradient_multishot(
    const std::vector<int32>& shot_indices,
    const std::vector<std::vector<std::vector<float64>>>& observed_data_all_shots
) {
    FwiGradientResult total_result;
    total_result.velocity_gradient = Array2D<float64>(nx_, nz_, 0.0);
    total_result.residual_norm = 0.0;
    total_result.num_time_steps = 0;

    for (size_t is = 0; is < shot_indices.size(); is++) {
        FwiGradientResult shot_result = compute_gradient(
            shot_indices[is],
            observed_data_all_shots[is]
        );

        for (int32 iz = 0; iz < nz_; iz++) {
            for (int32 ix = 0; ix < nx_; ix++) {
                total_result.velocity_gradient(ix, iz) += shot_result.velocity_gradient(ix, iz);
            }
        }

        total_result.residual_norm += shot_result.residual_norm;
        total_result.num_time_steps = shot_result.num_time_steps;
    }

    last_gradient_ = total_result.velocity_gradient;
    last_residual_ = total_result.residual_norm;

    return total_result;
}

}
