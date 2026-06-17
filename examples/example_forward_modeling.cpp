#include "fwi/fwi_engine.h"
#include <iostream>
#include <filesystem>
#include <chrono>

namespace fs = std::filesystem;

int main() {
    std::cout << "=========================================" << std::endl;
    std::cout << "  Seismic FWI Engine - Forward Modeling  " << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << std::endl;

    int nx = 200;
    int nz = 150;
    double dx = 10.0;
    double dz = 10.0;
    int num_sources = 3;
    int num_receivers = 100;
    double peak_freq = 25.0;
    int num_time_steps = 1500;
    int snapshot_interval = 50;
    int pml_width = 20;

    std::cout << "Model Parameters:" << std::endl;
    std::cout << "  Grid: " << nx << " x " << nz << " (dx=" << dx << "m, dz=" << dz << "m)" << std::endl;
    std::cout << "  Domain: " << (nx - 1) * dx << "m x " << (nz - 1) * dz << "m" << std::endl;
    std::cout << "  Sources: " << num_sources << ", Receivers: " << num_receivers << std::endl;
    std::cout << "  Peak frequency: " << peak_freq << " Hz" << std::endl;
    std::cout << "  Time steps: " << num_time_steps << std::endl;
    std::cout << "  PML width: " << pml_width << std::endl;
    std::cout << std::endl;

    fwi::FwiEngine engine;
    engine.set_compute_device(fwi::ComputeDevice::CPU);
    engine.set_vtk_format(fwi::VtkFormat::LEGACY_ASCII);

    std::cout << "Creating velocity model..." << std::endl;
    engine.create_synthetic_velocity_model(
        nx, nz, dx, dz,
        1500.0, 2500.0, 3500.0
    );

    const auto& model = engine.velocity_model();
    std::cout << "  Velocity range: " << model.get_min_velocity() << " - "
              << model.get_max_velocity() << " m/s" << std::endl;

    std::cout << "Creating acquisition geometry..." << std::endl;
    double source_x_start = pml_width * dx + 100.0;
    double source_x_end = (nx - 1 - pml_width) * dx - 100.0;
    double receiver_x_start = pml_width * dx + 50.0;
    double receiver_x_end = (nx - 1 - pml_width) * dx - 50.0;

    engine.create_survey(
        source_x_start, source_x_end, num_sources,
        receiver_x_start, receiver_x_end, num_receivers,
        50.0
    );

    std::cout << "  Shots: " << engine.geometry().num_shots() << std::endl;
    std::cout << "  Receivers: " << engine.geometry().num_receivers() << std::endl;

    std::cout << "Setting up wavelet..." << std::endl;
    double stable_dt = model.compute_stable_dt(0.5);
    std::cout << "  Stable dt (CFL=0.5): " << stable_dt * 1000 << " ms" << std::endl;

    fwi::Wavelet wavelet = fwi::Wavelet::ricker(stable_dt, peak_freq, stable_dt * 100);
    wavelet.normalize(1.0);
    engine.set_wavelet(wavelet);

    std::cout << "  Wavelet: Ricker, " << wavelet.num_samples() << " samples" << std::endl;

    fwi::FDSimulationParams params;
    params.dt = stable_dt;
    params.num_time_steps = num_time_steps;
    params.pml_width = pml_width;
    params.snapshot_interval = snapshot_interval;
    params.record_interval = 1;
    params.use_pml = true;
    params.record_receivers = true;
    engine.set_simulation_params(params);

    std::cout << "  Total simulation time: " << params.dt * params.num_time_steps << " s" << std::endl;
    std::cout << std::endl;

    std::string output_dir = "./output";
    fs::create_directories(output_dir);

    std::cout << "Exporting velocity model..." << std::endl;
    std::string vel_file = output_dir + "/velocity_model.vtk";
    if (engine.export_velocity_model(vel_file)) {
        std::cout << "  Saved to: " << vel_file << std::endl;
    }

    std::cout << std::endl;
    for (int shot_idx = 0; shot_idx < num_sources; shot_idx++) {
        std::cout << "=========================================" << std::endl;
        std::cout << "Processing shot " << shot_idx + 1 << "/" << num_sources << std::endl;
        std::cout << "=========================================" << std::endl;

        const auto& shot = engine.geometry().shot(shot_idx);
        std::cout << "  Shot location: (" << shot.x() << "m, " << shot.z() << "m)" << std::endl;

        auto start_time = std::chrono::high_resolution_clock::now();

        if (engine.run_forward_modeling(shot_idx)) {
            auto end_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);

            std::cout << "  Simulation completed in " << duration.count() << " seconds" << std::endl;
            std::cout << "  Snapshots recorded: " << engine.snapshots().size() << std::endl;

            std::string shot_dir = output_dir + "/shot_" + std::to_string(shot_idx);
            fs::create_directories(shot_dir);

            std::cout << "  Exporting snapshots..." << std::endl;
            if (engine.export_all_snapshots(shot_dir, "snapshot")) {
                std::cout << "  Snapshots exported to: " << shot_dir << std::endl;
            }

            std::string rec_file = shot_dir + "/receiver_data.txt";
            if (engine.export_receiver_data(rec_file)) {
                std::cout << "  Receiver data exported to: " << rec_file << std::endl;
            }

            if (!engine.snapshots().empty()) {
                const auto& last_snap = engine.snapshots().back();
                double max_p = 0.0;
                for (int64_t i = 0; i < last_snap.pressure.size(); i++) {
                    max_p = std::max(max_p, std::abs(last_snap.pressure[i]));
                }
                std::cout << "  Max pressure amplitude: " << max_p << std::endl;
            }

            if (!engine.receiver_data().empty()) {
                const auto& rec0 = engine.receiver_data()[0];
                double max_amp = 0.0;
                for (double v : rec0) {
                    max_amp = std::max(max_amp, std::abs(v));
                }
                std::cout << "  First receiver max amplitude: " << max_amp << std::endl;
                std::cout << "  Time samples per receiver: " << rec0.size() << std::endl;
            }
        } else {
            std::cout << "  ERROR: Forward modeling failed for shot " << shot_idx << std::endl;
            return 1;
        }

        std::cout << std::endl;
    }

    std::cout << "=========================================" << std::endl;
    std::cout << "  Forward modeling completed successfully!" << std::endl;
    std::cout << "  Output directory: " << fs::absolute(output_dir) << std::endl;
    std::cout << "=========================================" << std::endl;

    return 0;
}
