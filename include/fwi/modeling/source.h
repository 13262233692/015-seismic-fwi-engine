#pragma once

#include "fwi/common/types.h"
#include "fwi/modeling/wavelet.h"
#include "fwi/geometry/shot.h"
#include <vector>

namespace fwi {

enum class SourceType {
    PRESSURE,
    PARTICLE_VELOCITY_X,
    PARTICLE_VELOCITY_Z,
    EXPLOSIVE
};

class Source {
public:
    Source();
    Source(const Shot& shot, const Wavelet& wavelet, SourceType type = SourceType::EXPLOSIVE);

    const Shot& shot() const;
    Shot& shot();

    const Wavelet& wavelet() const;
    Wavelet& wavelet();

    SourceType type() const;
    void set_type(SourceType type);

    Point3D location() const;
    float64 x() const;
    float64 y() const;
    float64 z() const;

    float64 at(int32 time_idx) const;
    float64 at_time(float64 t) const;

    void add_to_grid(Array2D<float64>& grid, const Grid2D& grid_params, int32 time_idx, float64 dt) const;

    void add_to_grid_distributed(Array2D<float64>& grid, const Grid2D& grid_params, int32 time_idx, float64 dt, float64 sigma = 1.0) const;

    void get_grid_indices(const Grid2D& grid_params, int32& ix, int32& iz) const;

private:
    Shot shot_;
    Wavelet wavelet_;
    SourceType type_;
};

}
