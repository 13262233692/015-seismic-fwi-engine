#include "fwi/io/vtk_writer.h"
#include "fwi/modeling/velocity_model.h"
#include "fwi/modeling/acoustic_2d_fd.h"
#include "fwi/modeling/wavelet.h"
#include "fwi/modeling/source.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <filesystem>
#include <vector>

namespace fs = std::filesystem;

int main() {
    std::cout << "=== Testing VTK Writer ===" << std::endl;

    int passed = 0;
    int total = 0;

    std::string output_dir = "./test_vtk_output";
    fs::create_directories(output_dir);

    std::cout << "\n1. Testing velocity model VTK export..." << std::endl;
    total++;
    try {
        int nx = 50;
        int nz = 40;
        double dx = 10.0;
        double dz = 10.0;

        fwi::VelocityModel2D model(nx, nz, dx, dz);
        std::vector<std::pair<double, double>> layers = {
            {0.0, 1500.0},
            {100.0, 2500.0},
            {250.0, 3500.0}
        };
        model.set_layered_velocity(layers);
        model.set_uniform_density(2.0);

        fwi::VtkWriter writer;
        writer.set_format(fwi::VtkFormat::LEGACY_ASCII);

        std::string filename = output_dir + "/velocity_model.vtk";
        if (writer.write_velocity_model(filename, model)) {
            std::cout << "  PASS: Velocity model written to " << filename << std::endl;
            passed++;

            std::ifstream file(filename);
            std::string line;
            int line_count = 0;
            while (std::getline(file, line) && line_count < 10) {
                std::cout << "    " << line << std::endl;
                line_count++;
            }
        } else {
            std::cout << "  FAIL: Cannot write velocity model" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n2. Testing snapshot VTK export..." << std::endl;
    total++;
    try {
        int nx = 60;
        int nz = 50;
        double dx = 10.0;
        double dz = 10.0;

        fwi::Grid2D grid(nx, nz, dx, dz);

        fwi::Snapshot2D snapshot;
        snapshot.time_step = 100;
        snapshot.time = 0.5;
        snapshot.pressure = fwi::Array2D<double>(nx, nz, 0.0);
        snapshot.velocity_x = fwi::Array2D<double>(nx, nz, 0.0);
        snapshot.velocity_z = fwi::Array2D<double>(nx, nz, 0.0);

        double cx = nx * dx / 2.0;
        double cz = nz * dz / 2.0;
        for (int ix = 0; ix < nx; ix++) {
            for (int iz = 0; iz < nz; iz++) {
                double x = ix * dx;
                double z = iz * dz;
                double r = std::sqrt((x - cx) * (x - cx) + (z - cz) * (z - cz));
                double sigma = 100.0;
                snapshot.pressure(ix, iz) = std::exp(-r * r / (2 * sigma * sigma));
                snapshot.velocity_x(ix, iz) = (x - cx) / sigma * snapshot.pressure(ix, iz);
                snapshot.velocity_z(ix, iz) = (z - cz) / sigma * snapshot.pressure(ix, iz);
            }
        }

        fwi::VtkWriter writer;

        std::string ascii_file = output_dir + "/snapshot_ascii.vtk";
        writer.set_format(fwi::VtkFormat::LEGACY_ASCII);
        if (writer.write_snapshot(ascii_file, snapshot, grid)) {
            std::cout << "  PASS: Legacy ASCII snapshot written to " << ascii_file << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Cannot write Legacy ASCII snapshot" << std::endl;
        }

        total++;
        std::string bin_file = output_dir + "/snapshot_binary.vtk";
        writer.set_format(fwi::VtkFormat::LEGACY_BINARY);
        if (writer.write_snapshot(bin_file, snapshot, grid)) {
            std::cout << "  PASS: Legacy Binary snapshot written to " << bin_file << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Cannot write Legacy Binary snapshot" << std::endl;
        }

        total++;
        std::string xml_file = output_dir + "/snapshot_xml.vti";
        writer.set_format(fwi::VtkFormat::XML_IMAGE_DATA);
        if (writer.write_snapshot(xml_file, snapshot, grid)) {
            std::cout << "  PASS: XML ImageData snapshot written to " << xml_file << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Cannot write XML ImageData snapshot" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n3. Testing 3D snapshot export..." << std::endl;
    total++;
    try {
        int nx = 50;
        int nz = 40;
        double dx = 10.0;
        double dz = 10.0;

        fwi::Grid2D grid(nx, nz, dx, dz);

        fwi::Snapshot2D snapshot;
        snapshot.time_step = 50;
        snapshot.time = 0.25;
        snapshot.pressure = fwi::Array2D<double>(nx, nz, 0.0);
        snapshot.velocity_x = fwi::Array2D<double>(nx, nz, 0.0);
        snapshot.velocity_z = fwi::Array2D<double>(nx, nz, 0.0);

        for (int ix = 0; ix < nx; ix++) {
            for (int iz = 0; iz < nz; iz++) {
                snapshot.pressure(ix, iz) = std::sin(0.05 * ix) * std::cos(0.05 * iz);
            }
        }

        fwi::VtkWriter writer;
        writer.set_format(fwi::VtkFormat::LEGACY_ASCII);

        std::string filename = output_dir + "/snapshot_3d.vtk";
        if (writer.write_3d_snapshot(filename, snapshot, grid, 20.0)) {
            std::cout << "  PASS: 3D snapshot written to " << filename << std::endl;
            passed++;

            std::ifstream file(filename);
            std::string line;
            int line_count = 0;
            while (std::getline(file, line) && line_count < 15) {
                if (line.find("DATASET") != std::string::npos ||
                    line.find("DIMENSIONS") != std::string::npos ||
                    line.find("POINTS") != std::string::npos ||
                    line.find("SCALARS") != std::string::npos) {
                    std::cout << "    " << line << std::endl;
                }
                line_count++;
            }
        } else {
            std::cout << "  FAIL: Cannot write 3D snapshot" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n4. Testing snapshot series export..." << std::endl;
    total++;
    try {
        int nx = 40;
        int nz = 30;
        double dx = 10.0;
        double dz = 10.0;

        fwi::Grid2D grid(nx, nz, dx, dz);
        std::vector<fwi::Snapshot2D> snapshots;

        for (int t = 0; t < 5; t++) {
            fwi::Snapshot2D snap;
            snap.time_step = t * 20;
            snap.time = t * 0.1;
            snap.pressure = fwi::Array2D<double>(nx, nz, 0.0);
            snap.velocity_x = fwi::Array2D<double>(nx, nz, 0.0);
            snap.velocity_z = fwi::Array2D<double>(nx, nz, 0.0);

            double cx = nx * dx / 2.0;
            double cz = nz * dz / 2.0;
            double r0 = 50.0 + t * 30.0;

            for (int ix = 0; ix < nx; ix++) {
                for (int iz = 0; iz < nz; iz++) {
                    double x = ix * dx;
                    double z = iz * dz;
                    double r = std::sqrt((x - cx) * (x - cx) + (z - cz) * (z - cz));
                    double sigma = 40.0;
                    snap.pressure(ix, iz) = std::exp(-(r - r0) * (r - r0) / (2 * sigma * sigma));
                }
            }
            snapshots.push_back(snap);
        }

        fwi::VtkWriter writer;
        writer.set_format(fwi::VtkFormat::LEGACY_ASCII);

        std::string series_dir = output_dir + "/series";
        if (writer.write_snapshot_series(series_dir, "wave", snapshots, grid)) {
            std::cout << "  PASS: Snapshot series written to " << series_dir << std::endl;

            int file_count = 0;
            for (const auto& entry : fs::directory_iterator(series_dir)) {
                if (entry.path().extension() == ".vtk") {
                    file_count++;
                    std::cout << "    " << entry.path().filename().string() << std::endl;
                }
            }

            if (file_count == 5) {
                std::cout << "  PASS: Correct number of files in series" << std::endl;
                passed++;
            } else {
                std::cout << "  FAIL: Incorrect number of files, got " << file_count << std::endl;
            }
        } else {
            std::cout << "  FAIL: Cannot write snapshot series" << std::endl;
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
