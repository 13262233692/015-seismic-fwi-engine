#include "fwi/modeling/velocity_model.h"
#include "fwi/common/utils.h"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace fwi {

VelocityModel2D::VelocityModel2D()
    : grid_(), velocity_(), density_() {}

VelocityModel2D::VelocityModel2D(const Grid2D& grid)
    : grid_(grid)
    , velocity_(grid.nx, grid.nz, 1500.0)
    , density_(grid.nx, grid.nz, 1000.0) {}

VelocityModel2D::VelocityModel2D(int32 nx, int32 nz, float64 dx, float64 dz, float64 ox, float64 oz)
    : grid_(nx, nz, dx, dz, ox, oz)
    , velocity_(nx, nz, 1500.0)
    , density_(nx, nz, 1000.0) {}

const Grid2D& VelocityModel2D::grid() const { return grid_; }
Grid2D& VelocityModel2D::grid() { return grid_; }

const Array2D<float64>& VelocityModel2D::velocity() const { return velocity_; }
Array2D<float64>& VelocityModel2D::velocity() { return velocity_; }

const Array2D<float64>& VelocityModel2D::density() const { return density_; }
Array2D<float64>& VelocityModel2D::density() { return density_; }

float64 VelocityModel2D::get_velocity(int32 ix, int32 iz) const {
    ix = clamp(ix, 0, grid_.nx - 1);
    iz = clamp(iz, 0, grid_.nz - 1);
    return velocity_(ix, iz);
}

float64 VelocityModel2D::get_density(int32 ix, int32 iz) const {
    ix = clamp(ix, 0, grid_.nx - 1);
    iz = clamp(iz, 0, grid_.nz - 1);
    return density_(ix, iz);
}

void VelocityModel2D::set_velocity(int32 ix, int32 iz, float64 value) {
    if (ix >= 0 && ix < grid_.nx && iz >= 0 && iz < grid_.nz) {
        velocity_(ix, iz) = value;
    }
}

void VelocityModel2D::set_density(int32 ix, int32 iz, float64 value) {
    if (ix >= 0 && ix < grid_.nx && iz >= 0 && iz < grid_.nz) {
        density_(ix, iz) = value;
    }
}

void VelocityModel2D::set_uniform_velocity(float64 velocity) {
    velocity_.fill(velocity);
}

void VelocityModel2D::set_uniform_density(float64 density) {
    density_.fill(density);
}

void VelocityModel2D::set_gradient_velocity(float64 top_velocity, float64 bottom_velocity) {
    for (int32 iz = 0; iz < grid_.nz; iz++) {
        float64 frac = static_cast<float64>(iz) / (grid_.nz - 1);
        float64 vel = top_velocity + frac * (bottom_velocity - top_velocity);
        for (int32 ix = 0; ix < grid_.nx; ix++) {
            velocity_(ix, iz) = vel;
        }
    }
}

void VelocityModel2D::set_layered_velocity(const std::vector<std::pair<float64, float64>>& layers) {
    for (const auto& layer : layers) {
        float64 depth = layer.first;
        float64 vel = layer.second;

        int32 iz_start = nearest_index(depth, grid_.oz, grid_.dz);
        iz_start = clamp(iz_start, 0, grid_.nz - 1);

        int32 iz_end = grid_.nz;

        for (int32 iz = iz_start; iz < iz_end; iz++) {
            for (int32 ix = 0; ix < grid_.nx; ix++) {
                velocity_(ix, iz) = vel;
            }
        }
    }
}

void VelocityModel2D::add_rectangular_anomaly(float64 x1, float64 z1, float64 x2, float64 z2, float64 delta_v) {
    int32 ix1 = nearest_index(std::min(x1, x2), grid_.ox, grid_.dx);
    int32 ix2 = nearest_index(std::max(x1, x2), grid_.ox, grid_.dx);
    int32 iz1 = nearest_index(std::min(z1, z2), grid_.oz, grid_.dz);
    int32 iz2 = nearest_index(std::max(z1, z2), grid_.oz, grid_.dz);

    ix1 = clamp(ix1, 0, grid_.nx - 1);
    ix2 = clamp(ix2, 0, grid_.nx - 1);
    iz1 = clamp(iz1, 0, grid_.nz - 1);
    iz2 = clamp(iz2, 0, grid_.nz - 1);

    for (int32 iz = iz1; iz <= iz2; iz++) {
        for (int32 ix = ix1; ix <= ix2; ix++) {
            velocity_(ix, iz) += delta_v;
        }
    }
}

float64 VelocityModel2D::get_max_velocity() const {
    float64 max_v = velocity_(0, 0);
    for (int64 i = 0; i < velocity_.size(); i++) {
        if (velocity_[i] > max_v) {
            max_v = velocity_[i];
        }
    }
    return max_v;
}

float64 VelocityModel2D::get_min_velocity() const {
    float64 min_v = velocity_(0, 0);
    for (int64 i = 0; i < velocity_.size(); i++) {
        if (velocity_[i] < min_v) {
            min_v = velocity_[i];
        }
    }
    return min_v;
}

float64 VelocityModel2D::compute_stable_dt(float64 cfl) const {
    float64 vmax = get_max_velocity();
    float64 hmin = std::min(grid_.dx, grid_.dz);
    return cfl * hmin / vmax;
}

int32 VelocityModel2D::compute_pml_size(float64 freq_max, int32 npml) const {
    float64 vmin = get_min_velocity();
    float64 lambda_min = vmin / freq_max;
    int32 min_pml = static_cast<int32>(std::ceil(2.0 * lambda_min / std::min(grid_.dx, grid_.dz)));
    return std::max(npml, min_pml);
}

}
