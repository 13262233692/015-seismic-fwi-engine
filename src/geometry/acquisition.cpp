#include "fwi/geometry/acquisition.h"
#include <set>
#include <cmath>
#include <iostream>

namespace fwi {

AcquisitionGeometry::AcquisitionGeometry() {}

const std::vector<Shot>& AcquisitionGeometry::shots() const { return shots_; }
std::vector<Shot>& AcquisitionGeometry::shots() { return shots_; }
int32 AcquisitionGeometry::num_shots() const { return static_cast<int32>(shots_.size()); }
const Shot& AcquisitionGeometry::shot(int32 idx) const { return shots_[idx]; }
Shot& AcquisitionGeometry::shot(int32 idx) { return shots_[idx]; }

const std::vector<Receiver>& AcquisitionGeometry::receivers() const { return receivers_; }
std::vector<Receiver>& AcquisitionGeometry::receivers() { return receivers_; }
int32 AcquisitionGeometry::num_receivers() const { return static_cast<int32>(receivers_.size()); }
const Receiver& AcquisitionGeometry::receiver(int32 idx) const { return receivers_[idx]; }
Receiver& AcquisitionGeometry::receiver(int32 idx) { return receivers_[idx]; }

void AcquisitionGeometry::add_shot(const Shot& shot) { shots_.push_back(shot); }
void AcquisitionGeometry::add_receiver(const Receiver& receiver) { receivers_.push_back(receiver); }

void AcquisitionGeometry::clear_shots() { shots_.clear(); }
void AcquisitionGeometry::clear_receivers() { receivers_.clear(); }
void AcquisitionGeometry::clear() {
    shots_.clear();
    receivers_.clear();
}

float64 AcquisitionGeometry::apply_scalar(int32 value, int16 scalar) const {
    if (scalar == 0) {
        return static_cast<float64>(value);
    } else if (scalar > 0) {
        return static_cast<float64>(value) * static_cast<float64>(scalar);
    } else {
        return static_cast<float64>(value) / static_cast<float64>(-scalar);
    }
}

bool AcquisitionGeometry::load_from_segy(const std::string& segy_file) {
    SegyParser parser;
    if (!parser.open(segy_file)) {
        return false;
    }

    if (!parser.parse_all_traces()) {
        return false;
    }

    clear();

    const auto& traces = parser.traces();
    if (traces.empty()) {
        return true;
    }

    std::set<int32> source_ids;
    std::set<int32> receiver_ids;

    for (const auto& trace : traces) {
        int32 source_id = trace.header.energy_source_point_number;
        int32 receiver_id = trace.header.trace_number_field_record;

        int16 coord_scalar = trace.header.scalar_coordinates;
        int16 elev_scalar = trace.header.scalar_elevations;

        float64 sx = apply_scalar(trace.header.source_x, coord_scalar);
        float64 sy = apply_scalar(trace.header.source_y, coord_scalar);
        float64 sz = apply_scalar(trace.header.source_depth_below_surface, elev_scalar);

        float64 gx = apply_scalar(trace.header.group_x, coord_scalar);
        float64 gy = apply_scalar(trace.header.group_y, coord_scalar);
        float64 gz = apply_scalar(trace.header.receiver_group_elevation, elev_scalar);

        if (source_ids.find(source_id) == source_ids.end()) {
            source_ids.insert(source_id);
            shots_.emplace_back(source_id, sx, sy, sz);
        }

        if (receiver_ids.find(receiver_id) == receiver_ids.end()) {
            receiver_ids.insert(receiver_id);
            Receiver recv(receiver_id, gx, gy, gz);
            receivers_.push_back(recv);
        }

        for (auto& recv : receivers_) {
            if (recv.id() == receiver_id && recv.data().empty()) {
                recv.set_data(trace.data);
                break;
            }
        }
    }

    return true;
}

void AcquisitionGeometry::create_simple_survey(
    float64 source_x_start, float64 source_x_end, int32 num_sources,
    float64 receiver_x_start, float64 receiver_x_end, int32 num_receivers,
    float64 depth, float64 y
) {
    clear();

    float64 source_dx = (num_sources > 1) ? (source_x_end - source_x_start) / (num_sources - 1) : 0.0;
    for (int32 i = 0; i < num_sources; i++) {
        float64 x = source_x_start + i * source_dx;
        shots_.emplace_back(i, x, y, depth);
    }

    float64 receiver_dx = (num_receivers > 1) ? (receiver_x_end - receiver_x_start) / (num_receivers - 1) : 0.0;
    for (int32 i = 0; i < num_receivers; i++) {
        float64 x = receiver_x_start + i * receiver_dx;
        receivers_.emplace_back(i, x, y, depth);
    }
}

}
