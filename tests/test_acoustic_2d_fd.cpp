#include "fwi/modeling/acoustic_2d_fd.h"
#include "fwi/modeling/velocity_model.h"
#include "fwi/modeling/wavelet.h"
#include "fwi/modeling/source.h"
#include "fwi/geometry/receiver.h"
#include <iostream>
#include <cmath>
#include <vector>

int main() {
    std::cout << "=== Testing Acoustic 2D FD Solver ===" << std::endl;

    int passed = 0;
    int total = 0;

    std::cout << "\n1. Testing solver creation and initialization..." << std::endl;
    total++;
    try {
        int nx = 100;
        int nz = 80;
        double dx = 10.0;
        double dz = 10.0;

        fwi::VelocityModel2D model(nx, nz, dx, dz);
        model.set_uniform_velocity(2500.0);
        model.set_uniform_density(2.0);

        fwi::FDSimulationParams params;
        params.dt = model.compute_stable_dt(0.5);
        params.num_time_steps = 200;
        params.pml_width = 10;
        params.snapshot_interval = 50;
        params.use_pml = true;
        params.record_receivers = true;

        std::cout << "  Simulation parameters:" << std::endl;
        std::cout << "    dt: " << params.dt * 1000 << " ms" << std::endl;
        std::cout << "    Time steps: " << params.num_time_steps << std::endl;
        std::cout << "    Total time: " << params.dt * params.num_time_steps << " s" << std::endl;
        std::cout << "    PML width: " << params.pml_width << std::endl;

        fwi::Acoustic2DFD solver(model, params);

        fwi::Wavelet wavelet = fwi::Wavelet::ricker(params.dt, 20.0);
        fwi::Shot shot(0, nx * dx / 2.0, 0.0, 50.0);
        fwi::Source source(shot, wavelet, fwi::SourceType::EXPLOSIVE);
        solver.add_source(source);

        std::vector<fwi::Receiver> receivers;
        for (int i = 0; i < 20; i++) {
            double rx = 200.0 + i * 30.0;
            receivers.emplace_back(i, rx, 0.0, 10.0);
        }
        solver.set_receivers(receivers);

        solver.initialize();

        std::cout << "  Solver initialized successfully" << std::endl;
        std::cout << "    Current pressure max: " << solver.current_pressure().nx() << "x"
                  << solver.current_pressure().nz() << std::endl;

        std::cout << "  PASS: Solver creation and initialization works" << std::endl;
        passed++;
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n2. Testing solver time stepping..." << std::endl;
    total++;
    try {
        int nx = 80;
        int nz = 60;
        double dx = 10.0;
        double dz = 10.0;

        fwi::VelocityModel2D model(nx, nz, dx, dz);
        model.set_uniform_velocity(2500.0);
        model.set_uniform_density(2.0);

        fwi::FDSimulationParams params;
        params.dt = model.compute_stable_dt(0.5);
        params.num_time_steps = 100;
        params.pml_width = 10;
        params.snapshot_interval = 20;
        params.use_pml = true;
        params.record_receivers = true;

        fwi::Acoustic2DFD solver(model, params);

        fwi::Wavelet wavelet = fwi::Wavelet::ricker(params.dt, 20.0, params.dt * 30);
        fwi::Shot shot(0, nx * dx / 2.0, 0.0, 30.0);
        fwi::Source source(shot, wavelet, fwi::SourceType::EXPLOSIVE);
        solver.add_source(source);

        std::vector<fwi::Receiver> receivers;
        for (int i = 0; i < 10; i++) {
            double rx = 100.0 + i * 50.0;
            receivers.emplace_back(i, rx, 0.0, 10.0);
        }
        solver.set_receivers(receivers);

        solver.initialize();

        int progress = 0;
        solver.set_step_callback([&progress](int step, double time) {
            int new_progress = (step * 100) / 100;
            if (new_progress >= progress + 10) {
                progress = new_progress;
                std::cout << "    Progress: " << progress << "% (step " << step << ", time "
                          << time * 1000 << " ms)" << std::endl;
            }
        });

        solver.run();

        std::cout << "  Simulation completed:" << std::endl;
        std::cout << "    Final time step: " << solver.current_time_step() << std::endl;
        std::cout << "    Final time: " << solver.current_time() * 1000 << " ms" << std::endl;
        std::cout << "    Number of snapshots: " << solver.snapshots().size() << std::endl;
        std::cout << "    Number of receiver data: " << solver.receiver_data().size() << std::endl;

        if (solver.current_time_step() == params.num_time_steps &&
            solver.snapshots().size() >= 4 &&
            solver.receiver_data().size() == 10) {
            std::cout << "  PASS: Solver time stepping works correctly" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Solver time stepping results incorrect" << std::endl;
        }

        total++;
        if (!solver.snapshots().empty()) {
            const auto& last_snap = solver.snapshots().back();
            std::cout << "  Last snapshot:" << std::endl;
            std::cout << "    Time step: " << last_snap.time_step << std::endl;
            std::cout << "    Time: " << last_snap.time * 1000 << " ms" << std::endl;
            std::cout << "    Pressure grid: " << last_snap.pressure.nx() << "x" << last_snap.pressure.nz() << std::endl;

            double max_p = 0.0;
            for (int64_t i = 0; i < last_snap.pressure.size(); i++) {
                max_p = std::max(max_p, std::abs(last_snap.pressure[i]));
            }
            std::cout << "    Max pressure amplitude: " << max_p << std::endl;

            if (max_p > 0.0) {
                std::cout << "  PASS: Snapshot contains valid data" << std::endl;
                passed++;
            } else {
                std::cout << "  FAIL: Snapshot has zero amplitude" << std::endl;
            }
        }

        total++;
        if (!solver.receiver_data().empty()) {
            const auto& rec0 = solver.receiver_data()[0];
            double max_amp = 0.0;
            for (double v : rec0) {
                max_amp = std::max(max_amp, std::abs(v));
            }
            std::cout << "  Receiver 0 max amplitude: " << max_amp << std::endl;
            std::cout << "  Receiver 0 sample count: " << rec0.size() << std::endl;

            if (rec0.size() == params.num_time_steps && max_amp > 0.0) {
                std::cout << "  PASS: Receiver recording works" << std::endl;
                passed++;
            } else {
                std::cout << "  FAIL: Receiver recording incorrect" << std::endl;
            }
        }

    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n3. Testing layered model simulation..." << std::endl;
    total++;
    try {
        int nx = 120;
        int nz = 80;
        double dx = 10.0;
        double dz = 10.0;

        fwi::VelocityModel2D model(nx, nz, dx, dz);
        std::vector<std::pair<double, double>> layers = {
            {0.0, 1500.0},
            {200.0, 2500.0},
            {400.0, 3500.0}
        };
        model.set_layered_velocity(layers);
        model.set_uniform_density(2.0);

        fwi::FDSimulationParams params;
        params.dt = model.compute_stable_dt(0.4);
        params.num_time_steps = 300;
        params.pml_width = 15;
        params.snapshot_interval = 100;
        params.use_pml = true;

        fwi::Acoustic2DFD solver(model, params);

        fwi::Wavelet wavelet = fwi::Wavelet::ricker(params.dt, 15.0, params.dt * 50);
        fwi::Shot shot(0, nx * dx / 2.0, 0.0, 50.0);
        fwi::Source source(shot, wavelet, fwi::SourceType::EXPLOSIVE);
        solver.add_source(source);

        solver.initialize();
        solver.run();

        std::cout << "  Layered model simulation completed:" << std::endl;
        std::cout << "    Time steps: " << solver.current_time_step() << std::endl;
        std::cout << "    Snapshots: " << solver.snapshots().size() << std::endl;

        if (solver.snapshots().size() >= 2) {
            std::cout << "  PASS: Layered model simulation works" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Layered model simulation failed" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n=== Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << "/" << total << std::endl;

    if (passed == total) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests failed!" << std::endl;
        return 1;
    }
}
