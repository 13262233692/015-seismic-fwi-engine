#include "fwi/fwi_engine.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <iostream>

namespace fwi {

FwiEngine::FwiEngine()
    : compute_device_(ComputeDevice::CPU)
    , velocity_model_()
    , geometry_()
    , sim_params_()
    , wavelet_()
    , solver_(std::make_unique<Acoustic2DFD>())
    , vtk_writer_()
    , vtk_format_(VtkFormat::LEGACY_ASCII)
{
    wavelet_ = Wavelet::ricker(0.001, 30.0);
}

FwiEngine::~FwiEngine() {}

void FwiEngine::set_compute_device(ComputeDevice device) {
    compute_device_ = device;
}

ComputeDevice FwiEngine::compute_device() const {
    return compute_device_;
}

bool FwiEngine::load_segy_data(const std::string& segy_file) {
    return geometry_.load_from_segy(segy_file);
}

void FwiEngine::create_synthetic_velocity_model(int32 nx, int32 nz, float64 dx, float64 dz,
                                                float64 water_velocity,
                                                float64 sediment_velocity,
                                                float64 basement_velocity) {
    velocity_model_ = VelocityModel2D(nx, nz, dx, dz);

    int32 water_bottom_iz = nz / 4;
    int32 basement_iz = nz / 2;

    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 ix = 0; ix < nx; ix++) {
            if (iz < water_bottom_iz) {
                velocity_model_.set_velocity(ix, iz, water_velocity);
                velocity_model_.set_density(ix, iz, 1000.0);
            } else if (iz < basement_iz) {
                float64 frac = static_cast<float64>(iz - water_bottom_iz) / (basement_iz - water_bottom_iz);
                float64 vel = water_velocity + frac * (sediment_velocity - water_velocity);
                velocity_model_.set_velocity(ix, iz, vel);
                float64 rho = 1000.0 + frac * 600.0;
                velocity_model_.set_density(ix, iz, rho);
            } else {
                float64 frac = static_cast<float64>(iz - basement_iz) / (nz - basement_iz);
                float64 vel = sediment_velocity + frac * (basement_velocity - sediment_velocity);
                velocity_model_.set_velocity(ix, iz, vel);
                float64 rho = 1600.0 + frac * 800.0;
                velocity_model_.set_density(ix, iz, rho);
            }
        }
    }

    float64 anomaly_center_x = nx * dx * 0.6;
    float64 anomaly_center_z = basement_iz * dz + (nz - basement_iz) * dz * 0.3;
    float64 anomaly_size_x = nx * dx * 0.15;
    float64 anomaly_size_z = nz * dz * 0.1;
    velocity_model_.add_rectangular_anomaly(
        anomaly_center_x - anomaly_size_x / 2,
        anomaly_center_z - anomaly_size_z / 2,
        anomaly_center_x + anomaly_size_x / 2,
        anomaly_center_z + anomaly_size_z / 2,
        500.0
    );
}

void FwiEngine::set_velocity_model(const VelocityModel2D& model) {
    velocity_model_ = model;
}

const VelocityModel2D& FwiEngine::velocity_model() const {
    return velocity_model_;
}

void FwiEngine::create_survey(float64 source_x_start, float64 source_x_end, int32 num_sources,
                              float64 receiver_x_start, float64 receiver_x_end, int32 num_receivers,
                              float64 depth) {
    const Grid2D& grid = velocity_model_.grid();
    geometry_.create_simple_survey(
        source_x_start, source_x_end, num_sources,
        receiver_x_start, receiver_x_end, num_receivers,
        depth, 0.0
    );
}

const AcquisitionGeometry& FwiEngine::geometry() const {
    return geometry_;
}

AcquisitionGeometry& FwiEngine::geometry() {
    return geometry_;
}

void FwiEngine::set_simulation_params(const FDSimulationParams& params) {
    sim_params_ = params;
}

const FDSimulationParams& FwiEngine::simulation_params() const {
    return sim_params_;
}

void FwiEngine::set_wavelet(const Wavelet& wavelet) {
    wavelet_ = wavelet;
}

const Wavelet& FwiEngine::wavelet() const {
    return wavelet_;
}

bool FwiEngine::run_forward_modeling(int32 shot_index) {
    if (geometry_.num_shots() == 0) {
        std::cerr << "Error: No shots defined in survey geometry." << std::endl;
        return false;
    }

    if (shot_index < 0 || shot_index >= geometry_.num_shots()) {
        std::cerr << "Error: Shot index out of range." << std::endl;
        return false;
    }

    if (sim_params_.dt <= 0.0) {
        sim_params_.dt = velocity_model_.compute_stable_dt(0.5);
        std::cout << "Auto-computed dt = " << sim_params_.dt << "s based on CFL=0.5" << std::endl;
    }

    solver_ = std::make_unique<Acoustic2DFD>(velocity_model_, sim_params_);

    const Shot& shot = geometry_.shot(shot_index);
    Source source(shot, wavelet_, SourceType::EXPLOSIVE);
    solver_->add_source(source);

    solver_->set_receivers(geometry_.receivers());

    int32 total_steps = sim_params_.num_time_steps;
    int32 last_percent = -1;

    solver_->set_step_callback([&](int32 step, float64 time) {
        int32 percent = static_cast<int32>(100.0 * step / total_steps);
        if (percent != last_percent && percent % 10 == 0) {
            std::cout << "\rProgress: " << percent << "% (" << step << "/" << total_steps
                      << " steps, t=" << std::fixed << std::setprecision(3) << time << "s)" << std::flush;
            last_percent = percent;
        }
    });

    std::cout << "Starting forward modeling for shot " << shot_index
              << " at (" << shot.x() << ", " << shot.z() << ")..." << std::endl;
    std::cout << "Grid: " << velocity_model_.grid().nx << "x" << velocity_model_.grid().nz
              << ", dx=" << velocity_model_.grid().dx << "m, dz=" << velocity_model_.grid().dz << "m" << std::endl;
    std::cout << "Time steps: " << total_steps << ", dt=" << sim_params_.dt << "s" << std::endl;
    std::cout << "Total recording time: " << total_steps * sim_params_.dt << "s" << std::endl;
    std::cout << "PML width: " << sim_params_.pml_width << " grid points" << std::endl;
    std::cout << "Snapshot interval: every " << sim_params_.snapshot_interval << " steps" << std::endl;

    solver_->run();

    std::cout << std::endl << "Forward modeling completed. Saved "
              << solver_->snapshots().size() << " snapshots." << std::endl;

    return true;
}

const std::vector<Snapshot2D>& FwiEngine::snapshots() const {
    return solver_->snapshots();
}

const std::vector<std::vector<float64>>& FwiEngine::receiver_data() const {
    return solver_->receiver_data();
}

const std::vector<float64>& FwiEngine::time_axis() const {
    return solver_->time_axis();
}

bool FwiEngine::export_velocity_model(const std::string& filename) {
    vtk_writer_.set_format(vtk_format_);
    return vtk_writer_.write_velocity_model(filename, velocity_model_);
}

bool FwiEngine::export_snapshot(int32 snapshot_index, const std::string& filename) {
    if (snapshot_index < 0 || snapshot_index >= static_cast<int32>(solver_->snapshots().size())) {
        return false;
    }

    vtk_writer_.set_format(vtk_format_);
    return vtk_writer_.write_snapshot(
        filename,
        solver_->snapshots()[snapshot_index],
        velocity_model_.grid()
    );
}

bool FwiEngine::export_all_snapshots(const std::string& output_dir, const std::string& prefix) {
    std::filesystem::create_directories(output_dir);
    vtk_writer_.set_format(vtk_format_);
    return vtk_writer_.write_snapshot_series(
        output_dir, prefix,
        solver_->snapshots(),
        velocity_model_.grid()
    );
}

bool FwiEngine::export_receiver_data(const std::string& filename) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    const auto& recv_data = solver_->receiver_data();
    const auto& time = solver_->time_axis();
    const auto& receivers = geometry_.receivers();

    file << "# Time";
    for (size_t i = 0; i < receivers.size(); i++) {
        file << ", Receiver_" << i << "_x=" << receivers[i].x();
    }
    file << "\n";

    for (size_t t = 0; t < time.size(); t++) {
        file << std::scientific << std::setprecision(6) << time[t];
        for (size_t r = 0; r < recv_data.size(); r++) {
            if (t < recv_data[r].size()) {
                file << ", " << std::scientific << std::setprecision(6) << recv_data[r][t];
            } else {
                file << ", 0.0";
            }
        }
        file << "\n";
    }

    file.close();
    return true;
}

void FwiEngine::set_vtk_format(VtkFormat format) {
    vtk_format_ = format;
}

const Acoustic2DFD& FwiEngine::solver() const {
    return *solver_;
}

Acoustic2DFD& FwiEngine::solver() {
    return *solver_;
}

void FwiEngine::set_gradient_mode(GradientComputationMode mode) {
    gradient_computer_.set_computation_mode(mode);
}

GradientComputationMode FwiEngine::gradient_mode() const {
    return gradient_computer_.computation_mode();
}

void FwiEngine::set_gradient_threads(int32 num_threads) {
    gradient_computer_.set_num_threads(num_threads);
}

int32 FwiEngine::gradient_threads() const {
    return gradient_computer_.num_threads();
}

FwiGradientResult FwiEngine::compute_gradient(
    int32 shot_index,
    const std::vector<std::vector<float64>>& observed_data
) {
    gradient_computer_.set_velocity_model(velocity_model_);
    gradient_computer_.set_params(sim_params_);
    gradient_computer_.set_wavelet(wavelet_);
    gradient_computer_.set_geometry(geometry_);

    return gradient_computer_.compute_gradient(shot_index, observed_data);
}

FwiGradientResult FwiEngine::compute_gradient_multishot(
    const std::vector<int32>& shot_indices,
    const std::vector<std::vector<std::vector<float64>>>& observed_data_all_shots
) {
    gradient_computer_.set_velocity_model(velocity_model_);
    gradient_computer_.set_params(sim_params_);
    gradient_computer_.set_wavelet(wavelet_);
    gradient_computer_.set_geometry(geometry_);

    return gradient_computer_.compute_gradient_multishot(shot_indices, observed_data_all_shots);
}

bool FwiEngine::update_velocity_model(const Array2D<float64>& gradient, float64 step_length) {
    int32 nx = velocity_model_.grid().nx;
    int32 nz = velocity_model_.grid().nz;

    if (gradient.nx() != nx || gradient.nz() != nz) {
        std::cerr << "Error: Gradient dimensions mismatch." << std::endl;
        return false;
    }

    double max_grad = 0.0;
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 ix = 0; ix < nx; ix++) {
            max_grad = (std::max)(max_grad, std::abs(gradient(ix, iz)));
        }
    }

    if (max_grad < 1e-15) {
        std::cerr << "Warning: Gradient is zero, no update applied." << std::endl;
        return true;
    }

    double scale = step_length / max_grad;

    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 ix = 0; ix < nx; ix++) {
            double old_vel = velocity_model_.get_velocity(ix, iz);
            double new_vel = old_vel - scale * gradient(ix, iz);
            new_vel = (std::clamp)(new_vel, 1000.0, 6000.0);
            velocity_model_.set_velocity(ix, iz, new_vel);
        }
    }

    return true;
}

bool FwiEngine::export_gradient(const std::string& filename) {
    const Array2D<float64>& grad = gradient_computer_.last_gradient();
    if (grad.nx() == 0 || grad.nz() == 0) {
        std::cerr << "Error: No gradient available to export." << std::endl;
        return false;
    }

    vtk_writer_.set_format(vtk_format_);

    Grid2D grad_grid(grad.nx(), grad.nz(),
                     velocity_model_.grid().dx,
                     velocity_model_.grid().dz,
                     velocity_model_.grid().ox,
                     velocity_model_.grid().oz);

    Snapshot2D grad_snap;
    grad_snap.time_step = 0;
    grad_snap.time = 0.0;
    grad_snap.pressure = grad;
    grad_snap.velocity_x = Array2D<float64>(grad.nx(), grad.nz(), 0.0);
    grad_snap.velocity_z = Array2D<float64>(grad.nx(), grad.nz(), 0.0);

    return vtk_writer_.write_snapshot(filename, grad_snap, grad_grid);
}

const FwiGradient& FwiEngine::gradient_computer() const {
    return gradient_computer_;
}

}
