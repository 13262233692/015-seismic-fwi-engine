#pragma once

#include "fwi/common/types.h"
#include <vector>
#include <functional>

namespace fwi {

class VelocityModel2D {
public:
    VelocityModel2D();
    VelocityModel2D(const Grid2D& grid);
    VelocityModel2D(int32 nx, int32 nz, float64 dx, float64 dz, float64 ox = 0.0, float64 oz = 0.0);

    const Grid2D& grid() const;
    Grid2D& grid();

    const Array2D<float64>& velocity() const;
    Array2D<float64>& velocity();

    const Array2D<float64>& density() const;
    Array2D<float64>& density();

    float64 get_velocity(int32 ix, int32 iz) const;
    float64 get_density(int32 ix, int32 iz) const;

    void set_velocity(int32 ix, int32 iz, float64 value);
    void set_density(int32 ix, int32 iz, float64 value);

    void set_uniform_velocity(float64 velocity);
    void set_uniform_density(float64 density);

    void set_gradient_velocity(float64 top_velocity, float64 bottom_velocity);
    void set_layered_velocity(const std::vector<std::pair<float64, float64>>& layers);

    void add_rectangular_anomaly(float64 x1, float64 z1, float64 x2, float64 z2, float64 delta_v);

    float64 get_max_velocity() const;
    float64 get_min_velocity() const;

    float64 compute_stable_dt(float64 cfl = 0.5) const;

    int32 compute_pml_size(float64 freq_max, int32 npml = 20) const;

private:
    Grid2D grid_;
    Array2D<float64> velocity_;
    Array2D<float64> density_;
};

}
