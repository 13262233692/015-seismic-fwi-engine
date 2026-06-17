#include "fwi/modeling/source.h"
#include "fwi/common/utils.h"
#include <cmath>

namespace fwi {

Source::Source()
    : shot_()
    , wavelet_()
    , type_(SourceType::EXPLOSIVE) {}

Source::Source(const Shot& shot, const Wavelet& wavelet, SourceType type)
    : shot_(shot)
    , wavelet_(wavelet)
    , type_(type) {}

const Shot& Source::shot() const { return shot_; }
Shot& Source::shot() { return shot_; }

const Wavelet& Source::wavelet() const { return wavelet_; }
Wavelet& Source::wavelet() { return wavelet_; }

SourceType Source::type() const { return type_; }
void Source::set_type(SourceType type) { type_ = type; }

Point3D Source::location() const { return shot_.location(); }
float64 Source::x() const { return shot_.x(); }
float64 Source::y() const { return shot_.y(); }
float64 Source::z() const { return shot_.z(); }

float64 Source::at(int32 time_idx) const {
    return wavelet_.at(time_idx);
}

float64 Source::at_time(float64 t) const {
    return wavelet_.at_time(t);
}

void Source::get_grid_indices(const Grid2D& grid_params, int32& ix, int32& iz) const {
    ix = nearest_index(shot_.x(), grid_params.ox, grid_params.dx);
    iz = nearest_index(shot_.z(), grid_params.oz, grid_params.dz);
    ix = clamp(ix, 0, grid_params.nx - 1);
    iz = clamp(iz, 0, grid_params.nz - 1);
}

void Source::add_to_grid(Array2D<float64>& grid, const Grid2D& grid_params, int32 time_idx, float64 dt) const {
    int32 ix, iz;
    get_grid_indices(grid_params, ix, iz);

    float64 amplitude = wavelet_.at(time_idx);
    grid(ix, iz) += amplitude;
}

void Source::add_to_grid_distributed(Array2D<float64>& grid, const Grid2D& grid_params, int32 time_idx, float64 dt, float64 sigma) const {
    int32 ix, iz;
    get_grid_indices(grid_params, ix, iz);

    float64 amplitude = wavelet_.at(time_idx);

    int32 radius = static_cast<int32>(std::ceil(3.0 * sigma));
    float64 two_sigma2 = 2.0 * sigma * sigma;

    for (int32 dx = -radius; dx <= radius; dx++) {
        for (int32 dz = -radius; dz <= radius; dz++) {
            int32 jx = ix + dx;
            int32 jz = iz + dz;

            if (jx >= 0 && jx < grid_params.nx && jz >= 0 && jz < grid_params.nz) {
                float64 dist2 = static_cast<float64>(dx * dx + dz * dz);
                float64 weight = std::exp(-dist2 / two_sigma2);
                grid(jx, jz) += amplitude * weight;
            }
        }
    }
}

}
