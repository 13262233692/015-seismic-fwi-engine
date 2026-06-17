#include "fwi/geometry/receiver.h"

namespace fwi {

Receiver::Receiver() : id_(0), location_() {}

Receiver::Receiver(int32 id, const Point3D& location)
    : id_(id), location_(location) {}

Receiver::Receiver(int32 id, float64 x, float64 y, float64 z)
    : id_(id), location_(x, y, z) {}

int32 Receiver::id() const { return id_; }
const Point3D& Receiver::location() const { return location_; }
float64 Receiver::x() const { return location_.x; }
float64 Receiver::y() const { return location_.y; }
float64 Receiver::z() const { return location_.z; }

void Receiver::set_id(int32 id) { id_ = id; }
void Receiver::set_location(const Point3D& loc) { location_ = loc; }
void Receiver::set_location(float64 x, float64 y, float64 z) {
    location_.x = x;
    location_.y = y;
    location_.z = z;
}

const std::vector<float64>& Receiver::data() const { return data_; }
std::vector<float64>& Receiver::data() { return data_; }
void Receiver::set_data(const std::vector<float64>& data) { data_ = data; }

}
