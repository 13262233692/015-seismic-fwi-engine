#include "fwi/modeling/velocity_model.h"
#include "fwi/io/vtk_writer.h"
#include <iostream>
#include <cmath>
#include <filesystem>
#include <iomanip>

namespace fs = std::filesystem;

void print_model_info(const fwi::VelocityModel2D& model, const std::string& name) {
    std::cout << "  " << name << ":" << std::endl;
    std::cout << "    Grid: " << model.grid().nx << " x " << model.grid().nz << std::endl;
    std::cout << "    Cell size: " << model.grid().dx << "m x " << model.grid().dz << "m" << std::endl;
    std::cout << "    Domain: " << (model.grid().nx - 1) * model.grid().dx << "m x "
              << (model.grid().nz - 1) * model.grid().dz << "m" << std::endl;
    std::cout << "    Velocity range: " << std::fixed << std::setprecision(1)
              << model.get_min_velocity() << " - " << model.get_max_velocity() << " m/s" << std::endl;

    int mid_x = model.grid().nx / 2;
    std::cout << "    Velocity profile at x=" << mid_x * model.grid().dx << "m:" << std::endl;
    for (int iz = 0; iz < model.grid().nz; iz += 10) {
        double v = model.get_velocity(mid_x, iz);
        double depth = iz * model.grid().dz;
        std::cout << "      z=" << std::setw(5) << depth << "m: "
                  << std::setw(7) << std::fixed << std::setprecision(1) << v << " m/s" << std::endl;
    }
    std::cout << std::endl;
}

int main() {
    std::cout << "=========================================" << std::endl;
    std::cout << "  Velocity Model Building Example        " << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << std::endl;

    int nx = 200;
    int nz = 150;
    double dx = 10.0;
    double dz = 10.0;

    std::cout << "Base grid: " << nx << " x " << nz << " (dx=" << dx << "m, dz=" << dz << "m)" << std::endl;
    std::cout << "Domain size: " << (nx - 1) * dx << "m x " << (nz - 1) * dz << "m" << std::endl;
    std::cout << std::endl;

    std::string output_dir = "./velocity_models";
    fs::create_directories(output_dir);

    fwi::VtkWriter vtk_writer;
    vtk_writer.set_format(fwi::VtkFormat::LEGACY_ASCII);

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "1. Uniform Velocity Model" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    {
        fwi::VelocityModel2D model(nx, nz, dx, dz);
        model.set_uniform_velocity(3000.0);
        model.set_uniform_density(2.5);
        print_model_info(model, "Uniform Model");

        double stable_dt = model.compute_stable_dt(0.5);
        int pml_size = model.compute_pml_size(30.0, 20);
        std::cout << "    Stable dt (CFL=0.5): " << stable_dt * 1000 << " ms" << std::endl;
        std::cout << "    Recommended PML size for 30Hz: " << pml_size << " points" << std::endl;

        std::string filename = output_dir + "/model_uniform.vtk";
        if (vtk_writer.write_velocity_model(filename, model)) {
            std::cout << "    Saved to: " << filename << std::endl;
        }
    }
    std::cout << std::endl;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "2. Gradient Velocity Model" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    {
        fwi::VelocityModel2D model(nx, nz, dx, dz);
        model.set_gradient_velocity(1500.0, 4000.0);
        model.set_uniform_density(2.2);
        print_model_info(model, "Gradient Model");

        double stable_dt = model.compute_stable_dt(0.5);
        std::cout << "    Stable dt (CFL=0.5): " << stable_dt * 1000 << " ms" << std::endl;

        std::string filename = output_dir + "/model_gradient.vtk";
        if (vtk_writer.write_velocity_model(filename, model)) {
            std::cout << "    Saved to: " << filename << std::endl;
        }
    }
    std::cout << std::endl;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "3. Layered Velocity Model (Marine)" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    {
        fwi::VelocityModel2D model(nx, nz, dx, dz);

        std::vector<std::pair<double, double>> layers = {
            {0.0, 1500.0},
            {100.0, 1800.0},
            {200.0, 2200.0},
            {350.0, 2800.0},
            {500.0, 3500.0},
            {800.0, 4200.0},
            {1000.0, 4500.0}
        };
        model.set_layered_velocity(layers);

        std::vector<double> densities = {1.0, 1.1, 1.5, 2.0, 2.4, 2.6, 2.7};
        for (int ix = 0; ix < nx; ix++) {
            for (int iz = 0; iz < nz; iz++) {
                double z = iz * dz;
                double density = 2.0;
                for (int i = (int)layers.size() - 1; i >= 0; i--) {
                    if (z >= layers[i].first) {
                        density = densities[i];
                        break;
                    }
                }
                model.set_density(ix, iz, density);
            }
        }

        print_model_info(model, "Layered Model");

        double stable_dt = model.compute_stable_dt(0.5);
        std::cout << "    Stable dt (CFL=0.5): " << stable_dt * 1000 << " ms" << std::endl;

        std::cout << "    Layer boundaries:" << std::endl;
        for (const auto& layer : layers) {
            std::cout << "      Depth " << layer.first << "m: " << layer.second << " m/s" << std::endl;
        }

        std::string filename = output_dir + "/model_layered.vtk";
        if (vtk_writer.write_velocity_model(filename, model)) {
            std::cout << "    Saved to: " << filename << std::endl;
        }
    }
    std::cout << std::endl;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "4. Model with Anomalies (Salt Body)" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    {
        fwi::VelocityModel2D model(nx, nz, dx, dz);

        std::vector<std::pair<double, double>> layers = {
            {0.0, 1500.0},
            {100.0, 2000.0},
            {250.0, 2500.0},
            {500.0, 3000.0},
            {1000.0, 3500.0}
        };
        model.set_layered_velocity(layers);
        model.set_uniform_density(2.2);

        model.add_rectangular_anomaly(
            600.0, 200.0,
            1200.0, 600.0,
            1500.0
        );

        model.add_rectangular_anomaly(
            300.0, 700.0,
            500.0, 900.0,
            -500.0
        );

        print_model_info(model, "Anomaly Model");

        double stable_dt = model.compute_stable_dt(0.4);
        std::cout << "    Stable dt (CFL=0.4): " << stable_dt * 1000 << " ms" << std::endl;
        int pml_size = model.compute_pml_size(25.0, 20);
        std::cout << "    Recommended PML size for 25Hz: " << pml_size << " points" << std::endl;

        std::cout << "    Anomalies added:" << std::endl;
        std::cout << "      Salt body: +1500 m/s at (600,200)-(1200,600)" << std::endl;
        std::cout << "      Low-velocity zone: -500 m/s at (300,700)-(500,900)" << std::endl;

        std::string filename = output_dir + "/model_anomalies.vtk";
        if (vtk_writer.write_velocity_model(filename, model)) {
            std::cout << "    Saved to: " << filename << std::endl;
        }
    }
    std::cout << std::endl;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "5. Custom Velocity Model (Complex)" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    {
        fwi::VelocityModel2D model(nx, nz, dx, dz);
        model.set_gradient_velocity(1500.0, 4000.0);
        model.set_uniform_density(2.0);

        double cx = nx * dx / 2.0;
        double cz = nz * dz / 2.0;
        for (int ix = 0; ix < nx; ix++) {
            for (int iz = 0; iz < nz; iz++) {
                double x = ix * dx;
                double z = iz * dz;

                double v = model.get_velocity(ix, iz);

                double r1 = std::sqrt((x - cx) * (x - cx) + (z - cz + 200) * (z - cz + 200));
                if (r1 < 200.0) {
                    double factor = 1.0 + 0.5 * std::cos(M_PI * r1 / 200.0);
                    v *= factor;
                }

                double r2 = std::sqrt((x - cx + 400) * (x - cx + 400) + (z - cz - 100) * (z - cz - 100));
                if (r2 < 150.0) {
                    double factor = 1.0 - 0.3 * std::cos(M_PI * r2 / 150.0);
                    v *= factor;
                }

                v += std::sin(0.02 * x) * std::cos(0.01 * z) * 50.0;

                model.set_velocity(ix, iz, v);
            }
        }

        print_model_info(model, "Complex Model");

        double stable_dt = model.compute_stable_dt(0.5);
        std::cout << "    Stable dt (CFL=0.5): " << stable_dt * 1000 << " ms" << std::endl;

        std::cout << "    Features:" << std::endl;
        std::cout << "      - Background gradient (1500-4000 m/s)" << std::endl;
        std::cout << "      - Gaussian high-velocity anomaly (+50%)" << std::endl;
        std::cout << "      - Gaussian low-velocity anomaly (-30%)" << std::endl;
        std::cout << "      - Small-scale sinusoidal perturbations" << std::endl;

        std::string filename = output_dir + "/model_complex.vtk";
        if (vtk_writer.write_velocity_model(filename, model)) {
            std::cout << "    Saved to: " << filename << std::endl;
        }
    }
    std::cout << std::endl;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Stability Analysis Summary" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

    fwi::VelocityModel2D test_model(100, 100, 10.0, 10.0);
    test_model.set_uniform_velocity(3000.0);

    std::cout << "  Reference model: uniform 3000 m/s, 10m grid" << std::endl;
    std::cout << std::endl;
    std::cout << std::setw(10) << "CFL"
              << std::setw(15) << "dt (ms)"
              << std::setw(20) << "Max freq (Hz)" << std::endl;
    std::cout << std::string(45, '-') << std::endl;

    for (double cfl : {0.25, 0.4, 0.5, 0.6, 0.7}) {
        double dt = test_model.compute_stable_dt(cfl);
        double max_freq = 1.0 / (6.0 * dt);
        std::cout << std::setw(10) << cfl
                  << std::setw(15) << std::fixed << std::setprecision(2) << dt * 1000
                  << std::setw(20) << std::fixed << std::setprecision(1) << max_freq << std::endl;
    }

    std::cout << std::endl;
    std::cout << "  PML size recommendations:" << std::endl;
    std::cout << std::setw(15) << "Freq (Hz)"
              << std::setw(15) << "PML points"
              << std::setw(15) << "PML width (m)" << std::endl;
    std::cout << std::string(45, '-') << std::endl;

    for (double freq : {10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 50.0}) {
        int pml = test_model.compute_pml_size(freq, 20);
        std::cout << std::setw(15) << std::fixed << std::setprecision(1) << freq
                  << std::setw(15) << pml
                  << std::setw(15) << pml * test_model.grid().dx << std::endl;
    }

    std::cout << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "  Velocity model generation complete!    " << std::endl;
    std::cout << "  Output directory: " << fs::absolute(output_dir) << std::endl;
    std::cout << "=========================================" << std::endl;

    return 0;
}
