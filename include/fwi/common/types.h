#pragma once

#include <cstdint>
#include <vector>
#include <string>
#include <memory>
#include <complex>
#include <cmath>
#include <algorithm>

namespace fwi {

using float32 = float;
using float64 = double;
using int16 = std::int16_t;
using int32 = std::int32_t;
using int64 = std::int64_t;
using uint16 = std::uint16_t;
using uint32 = std::uint32_t;
using uint64 = std::uint64_t;

struct Point2D {
    float64 x;
    float64 z;

    Point2D() : x(0.0), z(0.0) {}
    Point2D(float64 x_, float64 z_) : x(x_), z(z_) {}
};

struct Point3D {
    float64 x;
    float64 y;
    float64 z;

    Point3D() : x(0.0), y(0.0), z(0.0) {}
    Point3D(float64 x_, float64 y_, float64 z_) : x(x_), y(y_), z(z_) {}
};

struct Grid2D {
    int32 nx;
    int32 nz;
    float64 dx;
    float64 dz;
    float64 ox;
    float64 oz;

    Grid2D() : nx(0), nz(0), dx(0.0), dz(0.0), ox(0.0), oz(0.0) {}
    Grid2D(int32 nx_, int32 nz_, float64 dx_, float64 dz_, float64 ox_ = 0.0, float64 oz_ = 0.0)
        : nx(nx_), nz(nz_), dx(dx_), dz(dz_), ox(ox_), oz(oz_) {}

    inline int64 index(int32 ix, int32 iz) const {
        return static_cast<int64>(iz) * nx + ix;
    }

    inline float64 coord_x(int32 ix) const {
        return ox + ix * dx;
    }

    inline float64 coord_z(int32 iz) const {
        return oz + iz * dz;
    }
};

template <typename T>
class Array2D {
public:
    Array2D() : nx_(0), nz_(0), data_() {}

    Array2D(int32 nx, int32 nz, T value = T())
        : nx_(nx), nz_(nz), data_(static_cast<size_t>(nx) * nz, value) {}

    Array2D(int32 nx, int32 nz, const std::vector<T>& data)
        : nx_(nx), nz_(nz), data_(data) {}

    inline int32 nx() const { return nx_; }
    inline int32 nz() const { return nz_; }
    inline int64 size() const { return static_cast<int64>(data_.size()); }

    inline T& operator()(int32 ix, int32 iz) {
        return data_[static_cast<size_t>(iz) * nx_ + ix];
    }

    inline const T& operator()(int32 ix, int32 iz) const {
        return data_[static_cast<size_t>(iz) * nx_ + ix];
    }

    inline T& operator[](int64 idx) { return data_[static_cast<size_t>(idx)]; }
    inline const T& operator[](int64 idx) const { return data_[static_cast<size_t>(idx)]; }

    inline T* data() { return data_.data(); }
    inline const T* data() const { return data_.data(); }

    inline std::vector<T>& vector() { return data_; }
    inline const std::vector<T>& vector() const { return data_; }

    void fill(T value) {
        std::fill(data_.begin(), data_.end(), value);
    }

private:
    int32 nx_;
    int32 nz_;
    std::vector<T> data_;
};

}
