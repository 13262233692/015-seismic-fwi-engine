#pragma once

#include "fwi/common/types.h"
#include "fwi/modeling/acoustic_2d_fd.h"
#include "fwi/modeling/velocity_model.h"
#include <string>
#include <vector>

namespace fwi {

enum class VtkFormat {
    LEGACY_ASCII,
    LEGACY_BINARY,
    XML_IMAGE_DATA
};

class VtkWriter {
public:
    VtkWriter();

    void set_format(VtkFormat format);
    VtkFormat format() const;

    bool write_snapshot(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid);

    bool write_velocity_model(const std::string& filename, const VelocityModel2D& model);

    bool write_snapshot_series(const std::string& output_dir, const std::string& prefix,
                               const std::vector<Snapshot2D>& snapshots, const Grid2D& grid);

    bool write_3d_snapshot(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid,
                           float64 y_thickness = 10.0);

private:
    VtkFormat format_;

    bool write_legacy_ascii(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid);
    bool write_legacy_binary(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid);
    bool write_xml_image_data(const std::string& filename, const Snapshot2D& snapshot, const Grid2D& grid);

    bool write_legacy_ascii_velocity(const std::string& filename, const VelocityModel2D& model);

    std::string format_number(float64 value) const;
};

}
