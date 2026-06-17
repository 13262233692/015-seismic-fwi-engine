#include "fwi/fwi_engine.h"
#include "fwi/modeling/fwi_gradient.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
#include <cstdint>

bool compare_gradients(const fwi::Array2D<double>& g1, const fwi::Array2D<double>& g2,
                       double& max_abs_diff, double& max_rel_diff) {
    if (g1.nx() != g2.nx() || g1.nz() != g2.nz()) {
        return false;
    }

    max_abs_diff = 0.0;
    max_rel_diff = 0.0;
    bool identical = true;

    for (int32_t iz = 0; iz < g1.nz(); iz++) {
        for (int32_t ix = 0; ix < g1.nx(); ix++) {
            double v1 = g1(ix, iz);
            double v2 = g2(ix, iz);
            double abs_diff = std::abs(v1 - v2);

            if (abs_diff > max_abs_diff) {
                max_abs_diff = abs_diff;
            }

            double denom = std::abs(v1) + std::abs(v2) + 1e-30;
            double rel_diff = abs_diff / denom;
            if (rel_diff > max_rel_diff) {
                max_rel_diff = rel_diff;
            }

            if (v1 != v2) {
                identical = false;
            }
        }
    }

    return identical;
}

uint64_t gradient_hash(const fwi::Array2D<double>& grad) {
    uint64_t hash = 0;
    const double* data = grad.data();
    int64_t n = static_cast<int64_t>(grad.nx()) * grad.nz();

    for (int64_t i = 0; i < n; i++) {
        uint64_t bits;
        std::memcpy(&bits, &data[i], sizeof(double));
        hash ^= bits + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
    }
    return hash;
}

int main() {
    std::cout << "================================================" << std::endl;
    std::cout << "  FWI Gradient Thread-Safety Test" << std::endl;
    std::cout << "================================================" << std::endl;
    std::cout << std::endl;

    const int32_t nx = 50;
    const int32_t nz = 40;
    const double dx = 20.0;
    const double dz = 20.0;
    const int32_t nt = 200;
    const double dt = 0.001;

    std::cout << "=== Test Configuration ===" << std::endl;
    std::cout << "  Grid: " << nx << " x " << nz << " (dx=" << dx << "m, dz=" << dz << "m)" << std::endl;
    std::cout << "  Time steps: " << nt << ", dt=" << dt << "s" << std::endl;
    std::cout << std::endl;

    fwi::FwiEngine engine;

    std::cout << "=== Creating Velocity Model ===" << std::endl;
    engine.create_synthetic_velocity_model(nx, nz, dx, dz, 1500.0, 2500.0, 3500.0);
    std::cout << "  Velocity range: " << engine.velocity_model().get_min_velocity()
              << " - " << engine.velocity_model().get_max_velocity() << " m/s" << std::endl;
    std::cout << std::endl;

    std::cout << "=== Creating Survey Geometry ===" << std::endl;
    engine.create_survey(
        200.0, 800.0, 3,
        100.0, 900.0, 20,
        20.0
    );
    std::cout << "  Shots: " << engine.geometry().num_shots() << std::endl;
    std::cout << "  Receivers: " << engine.geometry().num_receivers() << std::endl;
    std::cout << std::endl;

    std::cout << "=== Setting Simulation Parameters ===" << std::endl;
    fwi::FDSimulationParams params;
    params.dt = dt;
    params.num_time_steps = nt;
    params.pml_width = 10;
    params.snapshot_interval = 10;
    params.use_pml = true;
    params.record_receivers = true;
    engine.set_simulation_params(params);
    std::cout << "  dt=" << dt << "s, nt=" << nt << ", PML=" << params.pml_width << std::endl;
    std::cout << std::endl;

    std::cout << "=== Setting Wavelet ===" << std::endl;
    fwi::Wavelet wavelet = fwi::Wavelet::ricker(dt, 30.0);
    engine.set_wavelet(wavelet);
    std::cout << "  Ricker wavelet, 30 Hz, " << wavelet.num_samples() << " samples" << std::endl;
    std::cout << std::endl;

    std::cout << "=== Generating Observed Data (Synthetic) ===" << std::endl;
    std::cout << "  Creating true model (with anomaly) for observed data..." << std::endl;

    fwi::FwiEngine true_engine;
    true_engine.create_synthetic_velocity_model(nx, nz, dx, dz, 1500.0, 2600.0, 3600.0);
    true_engine.create_survey(
        200.0, 800.0, 3,
        100.0, 900.0, 20,
        20.0
    );
    true_engine.set_simulation_params(params);
    true_engine.set_wavelet(wavelet);
    true_engine.run_forward_modeling(0);
    auto observed_data = true_engine.receiver_data();
    std::cout << "  Generated " << observed_data.size() << " receiver traces, "
              << observed_data[0].size() << " samples each" << std::endl;
    std::cout << "  True model velocity range: "
              << true_engine.velocity_model().get_min_velocity() << " - "
              << true_engine.velocity_model().get_max_velocity() << " m/s" << std::endl;
    std::cout << "  Initial model velocity range: "
              << engine.velocity_model().get_min_velocity() << " - "
              << engine.velocity_model().get_max_velocity() << " m/s" << std::endl;
    std::cout << std::endl;

    const int32_t num_runs = 5;
    std::vector<std::vector<fwi::Array2D<double>>> results(4, std::vector<fwi::Array2D<double>>(num_runs));
    std::vector<std::string> mode_names = {"SERIAL", "OPENMP_UNSAFE", "OPENMP_ATOMIC", "OPENMP_REDUCTION"};

    std::cout << "=== Running Gradient Computations ===" << std::endl;
    std::cout << std::endl;

    for (int32_t mode = 0; mode < 4; mode++) {
        fwi::GradientComputationMode comp_mode;
        if (mode == 0) comp_mode = fwi::GradientComputationMode::SERIAL;
        else if (mode == 1) comp_mode = fwi::GradientComputationMode::OPENMP_UNSAFE;
        else if (mode == 2) comp_mode = fwi::GradientComputationMode::OPENMP_ATOMIC;
        else comp_mode = fwi::GradientComputationMode::OPENMP_REDUCTION;

        engine.set_gradient_mode(comp_mode);
        engine.set_gradient_threads(4);

        std::cout << "--- Mode: " << mode_names[mode] << " ---" << std::endl;

        for (int32_t run = 0; run < num_runs; run++) {
            fwi::FwiGradientResult result = engine.compute_gradient(0, observed_data);
            results[mode][run] = result.velocity_gradient;

            double grad_norm = 0.0;
            for (int32_t iz = 0; iz < nz; iz++) {
                for (int32_t ix = 0; ix < nx; ix++) {
                    double g = result.velocity_gradient(ix, iz);
                    grad_norm += g * g;
                }
            }
            grad_norm = std::sqrt(grad_norm);

            uint64_t hash = gradient_hash(result.velocity_gradient);

            std::cout << "  Run " << run + 1 << "/" << num_runs
                      << ": ||grad||=" << std::scientific << std::setprecision(6) << grad_norm
                      << ", hash=0x" << std::hex << hash << std::dec
                      << ", residual=" << result.residual_norm << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "=== Reproducibility Analysis ===" << std::endl;
    std::cout << std::endl;

    int32_t passed = 0;
    int32_t total = 0;

    for (int32_t mode = 0; mode < 4; mode++) {
        total++;
        std::cout << "--- " << mode_names[mode] << " ---" << std::endl;

        bool all_identical = true;
        double max_abs_diff = 0.0;
        double max_rel_diff = 0.0;

        for (int32_t run = 1; run < num_runs; run++) {
            double abs_diff, rel_diff;
            bool identical = compare_gradients(
                results[mode][0], results[mode][run],
                abs_diff, rel_diff
            );

            if (abs_diff > max_abs_diff) max_abs_diff = abs_diff;
            if (rel_diff > max_rel_diff) max_rel_diff = rel_diff;

            if (!identical) {
                all_identical = false;
            }
        }

        if (mode == 0 || mode == 3) {
            if (all_identical) {
                passed++;
                std::cout << "  ALL " << num_runs << " RUNS ARE BIT-IDENTICAL [OK]" << std::endl;
            } else {
                std::cout << "  RUN 0 vs 1: max_abs_diff=" << std::scientific << max_abs_diff
                          << ", max_rel_diff=" << std::scientific << max_rel_diff << std::endl;
                std::cout << "  FAIL: Runs are NOT reproducible [FAIL]" << std::endl;
            }
        } else {
            if (!all_identical) {
                std::cout << "  (expected) Runs are NOT bit-identical" << std::endl;
                std::cout << "  Max abs diff: " << std::scientific << max_abs_diff << std::endl;
                std::cout << "  Max rel diff: " << std::scientific << max_rel_diff << std::endl;
                if (mode == 1) {
                    if (max_rel_diff > 1e-3) {
                        passed++;
                        std::cout << "  (expected) Severe data race detected [OK]" << std::endl;
                    } else {
                        std::cout << "  Data race not severe enough [FAIL]" << std::endl;
                    }
                } else {
                    if (max_rel_diff < 1e-10) {
                        passed++;
                        std::cout << "  (expected) Only floating-point rounding differences [OK]" << std::endl;
                    } else {
                        std::cout << "  Unexpectedly large differences [FAIL]" << std::endl;
                    }
                }
            } else {
                std::cout << "  Unexpected: All runs are bit-identical [FAIL]" << std::endl;
            }
        }

        if (mode > 0) {
            double abs_diff, rel_diff;
            compare_gradients(
                results[0][0], results[mode][0],
                abs_diff, rel_diff
            );
            total++;
            std::cout << "  vs SERIAL: max_abs_diff=" << std::scientific << abs_diff
                      << ", max_rel_diff=" << std::scientific << rel_diff;
            if (mode == 1) {
                if (rel_diff > 1e-3) {
                    passed++;
                    std::cout << " (data race corruption) [OK]" << std::endl;
                } else {
                    std::cout << " [FAIL]" << std::endl;
                }
            } else {
                if (abs_diff < 1e-12) {
                    passed++;
                    std::cout << " [OK]" << std::endl;
                } else {
                    std::cout << " [FAIL]" << std::endl;
                }
            }
        }
        std::cout << std::endl;
    }

    std::cout << "=== Summary ===" << std::endl;
    std::cout << "  Passed: " << passed << "/" << total << std::endl;
    std::cout << std::endl;

    std::cout << "=== Exporting Results ===" << std::endl;
    engine.export_gradient("gradient_serial.vtk");
    std::cout << "  Saved gradient to gradient_serial.vtk" << std::endl;
    std::cout << std::endl;

    if (passed == total) {
        std::cout << "================================================" << std::endl;
        std::cout << "  ALL TESTS PASSED [OK]" << std::endl;
        std::cout << "================================================" << std::endl;
        return 0;
    } else {
        std::cout << "================================================" << std::endl;
        std::cout << "  SOME TESTS FAILED [FAIL]" << std::endl;
        std::cout << "================================================" << std::endl;
        return 1;
    }
}
