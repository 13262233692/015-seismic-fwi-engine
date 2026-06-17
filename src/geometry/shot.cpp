#include "fwi/geometry/shot.h"

namespace fwi {

Shot::Shot() : id_(0), location_() {}

Shot::Shot(int32 id, const Point3D& location)
    : id_(id), location_(location) {}

Shot::Shot(int32 id, float64 x, float64 y, float64 z)
    : id_(id), location_(x, y, z) {}

int32 Shot::id() const { return id_; }
const Point3D& Shot::location() const { return location_; }
float64 Shot::x() const { return location_.x; }
float64 Shot::y() const { return location_.y; }
float64 Shot::z() const { return location_.z; }

void Shot::set_id(int32 id) { id_ = id; }
void Shot::set_location(const Point3D& loc) { location_ = loc; }
void Shot::set_location(float64 x, float64 y, float64 z) {
    location_.x = x;
    location_.y = y;
    location_.z = z;
}

}
