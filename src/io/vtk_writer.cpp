#include "fwi/io/vtk_writer.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <filesystem>
#include <stdexcept>

namespace fwi {

VtkWriter::VtkWriter() : format_(VtkFormat::LEGACY_ASCII) {}

void VtkWriter::set_format(VtkFormat format) { format_ = format; }
VtkFormat VtkWriter::format() const { return format_; }

std::string VtkWriter::format_number(float64 value) const {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(6) << value;
    return oss.str();
}

bool VtkWriter::write_snapshot(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid) {
    switch (format_) {
        case VtkFormat::LEGACY_ASCII:
            return write_legacy_ascii(filename, snapshot, grid);
        case VtkFormat::LEGACY_BINARY:
            return write_legacy_binary(filename, snapshot, grid);
        case VtkFormat::XML_IMAGE_DATA:
            return write_xml_image_data(filename, snapshot, grid);
        default:
            return false;
    }
}

bool VtkWriter::write_velocity_model(const std::string& filename, const VelocityModel2D& model) {
    return write_legacy_ascii_velocity(filename, model);
}

bool VtkWriter::write_legacy_ascii(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    int32 nx = grid.nx;
    int32 nz = grid.nz;
    int32 ny = 1;
    float64 dx = grid.dx;
    float64 dz = grid.dz;
    float64 dy = std::min(dx, dz);
    float64 ox = grid.ox;
    float64 oz = grid.oz;
    float64 oy = 0.0;

    file << "# vtk DataFile Version 3.0\n";
    file << "Seismic Wavefield Snapshot - Time: " << snapshot.time << "s, Step: " << snapshot.time_step << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    file << "ORIGIN " << format_number(ox) << " " << format_number(oy) << " " << format_number(oz) << "\n";
    file << "SPACING " << format_number(dx) << " " << format_number(dy) << " " << format_number(dz) << "\n";
    file << "POINT_DATA " << nx * ny * nz << "\n";

    file << "SCALARS Pressure float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                file << format_number(snapshot.pressure(ix, iz)) << " ";
            }
            file << "\n";
        }
    }

    file << "SCALARS VelocityX float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                file << format_number(snapshot.velocity_x(ix, iz)) << " ";
            }
            file << "\n";
        }
    }

    file << "SCALARS VelocityZ float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                file << format_number(snapshot.velocity_z(ix, iz)) << " ";
            }
            file << "\n";
        }
    }

    float64 vmax = 0.0;
    for (int64 i = 0; i < snapshot.pressure.size(); i++) {
        vmax = std::max(vmax, std::abs(snapshot.pressure[i]));
    }
    file << "FIELD FieldData 2\n";
    file << "TIME 1 1 float\n";
    file << format_number(snapshot.time) << "\n";
    file << "TIMESTEP 1 1 int\n";
    file << snapshot.time_step << "\n";
    file << "MAX_PRESSURE 1 1 float\n";
    file << format_number(vmax) << "\n";

    file.close();
    return true;
}

bool VtkWriter::write_legacy_binary(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid) {
    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        return false;
    }

    int32 nx = grid.nx;
    int32 nz = grid.nz;
    int32 ny = 1;
    float64 dx = grid.dx;
    float64 dz = grid.dz;
    float64 dy = std::min(dx, dz);
    float64 ox = grid.ox;
    float64 oz = grid.oz;
    float64 oy = 0.0;

    std::ostringstream header;
    header << "# vtk DataFile Version 3.0\n";
    header << "Seismic Wavefield Snapshot - Time: " << snapshot.time << "s, Step: " << snapshot.time_step << "\n";
    header << "BINARY\n";
    header << "DATASET STRUCTURED_POINTS\n";
    header << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    header << "ORIGIN " << ox << " " << oy << " " << oz << "\n";
    header << "SPACING " << dx << " " << dy << " " << dz << "\n";
    header << "POINT_DATA " << nx * ny * nz << "\n";

    file << header.str();

    file << "SCALARS Pressure float 1\n";
    file << "LOOKUP_TABLE default\n";
    std::vector<float32> p_data(snapshot.pressure.size());
    for (int64 i = 0; i < snapshot.pressure.size(); i++) {
        p_data[i] = static_cast<float32>(snapshot.pressure[i]);
    }
    file.write(reinterpret_cast<char*>(p_data.data()), p_data.size() * sizeof(float32));

    file << "\nSCALARS VelocityX float 1\n";
    file << "LOOKUP_TABLE default\n";
    std::vector<float32> vx_data(snapshot.velocity_x.size());
    for (int64 i = 0; i < snapshot.velocity_x.size(); i++) {
        vx_data[i] = static_cast<float32>(snapshot.velocity_x[i]);
    }
    file.write(reinterpret_cast<char*>(vx_data.data()), vx_data.size() * sizeof(float32));

    file << "\nSCALARS VelocityZ float 1\n";
    file << "LOOKUP_TABLE default\n";
    std::vector<float32> vz_data(snapshot.velocity_z.size());
    for (int64 i = 0; i < snapshot.velocity_z.size(); i++) {
        vz_data[i] = static_cast<float32>(snapshot.velocity_z[i]);
    }
    file.write(reinterpret_cast<char*>(vz_data.data()), vz_data.size() * sizeof(float32));

    file.close();
    return true;
}

bool VtkWriter::write_xml_image_data(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    int32 nx = grid.nx;
    int32 nz = grid.nz;
    int32 ny = 1;
    float64 dx = grid.dx;
    float64 dz = grid.dz;
    float64 dy = std::min(dx, dz);
    float64 ox = grid.ox;
    float64 oz = grid.oz;
    float64 oy = 0.0;

    file << "<?xml version=\"1.0\"?>\n";
    file << "<VTKFile type=\"ImageData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
    file << "  <ImageData WholeExtent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\"\n";
    file << "             Origin=\"" << ox << " " << oy << " " << oz << "\"\n";
    file << "             Spacing=\"" << dx << " " << dy << " " << dz << "\">\n";
    file << "    <Piece Extent=\"0 " << nx - 1 << " 0 " << ny - 1 << " 0 " << nz - 1 << "\">\n";
    file << "      <PointData Scalars=\"Pressure\">\n";

    file << "        <DataArray type=\"Float32\" Name=\"Pressure\" format=\"ascii\">\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                file << " " << snapshot.pressure(ix, iz);
            }
        }
    }
    file << "\n        </DataArray>\n";

    file << "        <DataArray type=\"Float32\" Name=\"VelocityX\" format=\"ascii\">\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                file << " " << snapshot.velocity_x(ix, iz);
            }
        }
    }
    file << "\n        </DataArray>\n";

    file << "        <DataArray type=\"Float32\" Name=\"VelocityZ\" format=\"ascii\">\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                file << " " << snapshot.velocity_z(ix, iz);
            }
        }
    }
    file << "\n        </DataArray>\n";

    file << "      </PointData>\n";
    file << "      <CellData>\n";
    file << "      </CellData>\n";
    file << "    </Piece>\n";
    file << "  </ImageData>\n";

    file << "  <FieldData>\n";
    file << "    <DataArray type=\"Float32\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">\n";
    file << "      " << snapshot.time << "\n";
    file << "    </DataArray>\n";
    file << "    <DataArray type=\"Int32\" Name=\"TIMESTEP\" NumberOfTuples=\"1\" format=\"ascii\">\n";
    file << "      " << snapshot.time_step << "\n";
    file << "    </DataArray>\n";
    file << "  </FieldData>\n";

    file << "</VTKFile>\n";

    file.close();
    return true;
}

bool VtkWriter::write_legacy_ascii_velocity(const std::string& filename, const VelocityModel2D& model) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    const Grid2D& grid = model.grid();
    int32 nx = grid.nx;
    int32 nz = grid.nz;
    int32 ny = 1;
    float64 dx = grid.dx;
    float64 dz = grid.dz;
    float64 dy = std::min(dx, dz);
    float64 ox = grid.ox;
    float64 oz = grid.oz;
    float64 oy = 0.0;

    file << "# vtk DataFile Version 3.0\n";
    file << "Velocity Model\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    file << "ORIGIN " << format_number(ox) << " " << format_number(oy) << " " << format_number(oz) << "\n";
    file << "SPACING " << format_number(dx) << " " << format_number(dy) << " " << format_number(dz) << "\n";
    file << "POINT_DATA " << nx * ny * nz << "\n";

    file << "SCALARS Velocity float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                file << format_number(model.get_velocity(ix, iz)) << " ";
            }
            file << "\n";
        }
    }

    file << "SCALARS Density float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                file << format_number(model.get_density(ix, iz)) << " ";
            }
            file << "\n";
        }
    }

    file.close();
    return true;
}

bool VtkWriter::write_snapshot_series(const std::string& output_dir, const std::string& prefix,
                                      const std::vector<Snapshot2D>& snapshots, const Grid2D& grid) {
    std::filesystem::create_directories(output_dir);

    for (size_t i = 0; i < snapshots.size(); i++) {
        std::ostringstream oss;
        oss << output_dir << "/" << prefix << "_" << std::setw(6) << std::setfill('0') << i << ".vtk";
        std::string filename = oss.str();

        if (!write_snapshot(filename, snapshots[i], grid)) {
            return false;
        }
    }

    return true;
}

bool VtkWriter::write_3d_snapshot(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid,
                                  float64 y_thickness) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        return false;
    }

    int32 nx = grid.nx;
    int32 nz = grid.nz;
    int32 ny = 3;
    float64 dx = grid.dx;
    float64 dz = grid.dz;
    float64 dy = y_thickness / (ny - 1);
    float64 ox = grid.ox;
    float64 oz = grid.oz;
    float64 oy = -y_thickness / 2.0;

    file << "# vtk DataFile Version 3.0\n";
    file << "3D Seismic Wavefield Snapshot - Time: " << snapshot.time << "s, Step: " << snapshot.time_step << "\n";
    file << "ASCII\n";
    file << "DATASET STRUCTURED_POINTS\n";
    file << "DIMENSIONS " << nx << " " << ny << " " << nz << "\n";
    file << "ORIGIN " << format_number(ox) << " " << format_number(oy) << " " << format_number(oz) << "\n";
    file << "SPACING " << format_number(dx) << " " << format_number(dy) << " " << format_number(dz) << "\n";
    file << "POINT_DATA " << nx * ny * nz << "\n";

    file << "SCALARS Pressure float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                float64 taper = 1.0;
                if (iy == 0 || iy == ny - 1) {
                    taper = 0.3;
                }
                file << format_number(snapshot.pressure(ix, iz) * taper) << " ";
            }
            file << "\n";
        }
    }

    file << "SCALARS VelocityX float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                float64 taper = 1.0;
                if (iy == 0 || iy == ny - 1) {
                    taper = 0.3;
                }
                file << format_number(snapshot.velocity_x(ix, iz) * taper) << " ";
            }
            file << "\n";
        }
    }

    file << "SCALARS VelocityZ float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int32 iz = 0; iz < nz; iz++) {
        for (int32 iy = 0; iy < ny; iy++) {
            for (int32 ix = 0; ix < nx; ix++) {
                float64 taper = 1.0;
                if (iy == 0 || iy == ny - 1) {
                    taper = 0.3;
                }
                file << format_number(snapshot.velocity_z(ix, iz) * taper) << " ";
            }
            file << "\n";
        }
    }

    file.close();
    return true;
}

}
