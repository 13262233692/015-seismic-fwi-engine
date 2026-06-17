#pragma once

#include "fwi/common/types.h"
#include <vector>

namespace fwi {

class Shot {
public:
    Shot();
    Shot(int32 id, const Point3D& location);
    Shot(int32 id, float64 x, float64 y, float64 z);

    int32 id() const;
    const Point3D& location() const;
    float64 x() const;
    float64 y() const;
    float64 z() const;

    void set_id(int32 id);
    void set_location(const Point3D& loc);
    void set_location(float64 x, float64 y, float64 z);

private:
    int32 id_;
    Point3D location_;
};

}
