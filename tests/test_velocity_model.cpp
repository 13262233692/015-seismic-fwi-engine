#include "fwi/modeling/velocity_model.h"
#include <iostream>
#include <cmath>
#include <vector>

int main() {
    std::cout << "=== Testing Velocity Model ===" << std::endl;

    int passed = 0;
    int total = 0;

    std::cout << "\n1. Testing basic velocity model creation..." << std::endl;
    total++;
    try {
        int nx = 100;
        int nz = 80;
        double dx = 10.0;
        double dz = 10.0;

        fwi::VelocityModel2D model(nx, nz, dx, dz);

        std::cout << "  Model grid:" << std::endl;
        std::cout << "    nx=" << model.grid().nx << ", nz=" << model.grid().nz << std::endl;
        std::cout << "    dx=" << model.grid().dx << "m, dz=" << model.grid().dz << "m" << std::endl;
        std::cout << "    ox=" << model.grid().ox << "m, oz=" << model.grid().oz << "m" << std::endl;

        if (model.grid().nx == nx && model.grid().nz == nz) {
            std::cout << "  PASS: Model dimensions correct" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Model dimensions incorrect" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n2. Testing uniform velocity model..." << std::endl;
    total++;
    try {
        int nx = 50;
        int nz = 50;
        fwi::VelocityModel2D model(nx, nz, 5.0, 5.0);

        double velocity = 2500.0;
        model.set_uniform_velocity(velocity);
        model.set_uniform_density(2.0);

        double max_v = model.get_max_velocity();
        double min_v = model.get_min_velocity();

        std::cout << "  Uniform velocity: " << velocity << " m/s" << std::endl;
        std::cout << "  Max velocity in model: " << max_v << " m/s" << std::endl;
        std::cout << "  Min velocity in model: " << min_v << " m/s" << std::endl;

        if (std::abs(max_v - velocity) < 1e-6 && std::abs(min_v - velocity) < 1e-6) {
            std::cout << "  PASS: Uniform velocity set correctly" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Uniform velocity not set correctly" << std::endl;
        }

        total++;
        double v_center = model.get_velocity(nx / 2, nz / 2);
        if (std::abs(v_center - velocity) < 1e-6) {
            std::cout << "  PASS: Velocity at center correct" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Velocity at center incorrect: " << v_center << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n3. Testing gradient velocity model..." << std::endl;
    total++;
    try {
        int nx = 100;
        int nz = 100;
        fwi::VelocityModel2D model(nx, nz, 5.0, 5.0);

        double top_v = 1500.0;
        double bottom_v = 3500.0;
        model.set_gradient_velocity(top_v, bottom_v);

        double v_top = model.get_velocity(50, 0);
        double v_bottom = model.get_velocity(50, nz - 1);
        double v_mid = model.get_velocity(50, nz / 2);

        std::cout << "  Velocity at top: " << v_top << " m/s" << std::endl;
        std::cout << "  Velocity at bottom: " << v_bottom << " m/s" << std::endl;
        std::cout << "  Velocity at middle: " << v_mid << " m/s" << std::endl;

        if (std::abs(v_top - top_v) < 1.0 && std::abs(v_bottom - bottom_v) < 1.0 && v_mid > top_v && v_mid < bottom_v) {
            std::cout << "  PASS: Gradient velocity model created correctly" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Gradient velocity model incorrect" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n4. Testing layered velocity model..." << std::endl;
    total++;
    try {
        int nx = 100;
        int nz = 100;
        fwi::VelocityModel2D model(nx, nz, 5.0, 5.0);

        std::vector<std::pair<double, double>> layers = {
            {0.0, 1500.0},
            {100.0, 2000.0},
            {200.0, 2500.0},
            {300.0, 3000.0}
        };
        model.set_layered_velocity(layers);

        double v_water = model.get_velocity(50, 5);
        double v_layer2 = model.get_velocity(50, 25);
        double v_layer3 = model.get_velocity(50, 50);

        std::cout << "  Velocity in water layer: " << v_water << " m/s" << std::endl;
        std::cout << "  Velocity in layer 2: " << v_layer2 << " m/s" << std::endl;
        std::cout << "  Velocity in layer 3: " << v_layer3 << " m/s" << std::endl;

        if (std::abs(v_water - 1500.0) < 10.0 && std::abs(v_layer2 - 2000.0) < 10.0) {
            std::cout << "  PASS: Layered velocity model created correctly" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Layered velocity model incorrect" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n5. Testing stability analysis..." << std::endl;
    total++;
    try {
        int nx = 200;
        int nz = 150;
        double dx = 10.0;
        double dz = 10.0;
        fwi::VelocityModel2D model(nx, nz, dx, dz);
        model.set_uniform_velocity(3000.0);

        double stable_dt = model.compute_stable_dt(0.5);
        std::cout << "  CFL=0.5 stable dt: " << stable_dt * 1000 << " ms" << std::endl;

        int pml_size = model.compute_pml_size(30.0, 20);
        std::cout << "  PML size for 30Hz: " << pml_size << " grid points" << std::endl;

        double min_expected_dt = 0.0;
        double max_expected_dt = 0.01;
        if (stable_dt > min_expected_dt && stable_dt < max_expected_dt && pml_size >= 10) {
            std::cout << "  PASS: Stability analysis results reasonable" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Stability analysis results out of expected range" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n6. Testing rectangular anomaly..." << std::endl;
    total++;
    try {
        int nx = 100;
        int nz = 80;
        fwi::VelocityModel2D model(nx, nz, 5.0, 5.0);
        model.set_uniform_velocity(2500.0);

        model.add_rectangular_anomaly(100.0, 100.0, 300.0, 250.0, 500.0);

        double v_outside = model.get_velocity(10, 10);
        double v_inside = model.get_velocity(40, 30);

        std::cout << "  Velocity outside anomaly: " << v_outside << " m/s" << std::endl;
        std::cout << "  Velocity inside anomaly: " << v_inside << " m/s" << std::endl;

        if (std::abs(v_outside - 2500.0) < 1.0 && v_inside > 2900.0) {
            std::cout << "  PASS: Rectangular anomaly added correctly" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Rectangular anomaly not added correctly" << std::endl;
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
