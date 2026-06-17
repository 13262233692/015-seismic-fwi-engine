#pragma once

#include "fwi/common/types.h"
#include "fwi/geometry/shot.h"
#include "fwi/geometry/receiver.h"
#include "fwi/segy/segy_parser.h"
#include <vector>
#include <memory>

namespace fwi {

class AcquisitionGeometry {
public:
    AcquisitionGeometry();

    const std::vector<Shot>& shots() const;
    std::vector<Shot>& shots();
    int32 num_shots() const;
    const Shot& shot(int32 idx) const;
    Shot& shot(int32 idx);

    const std::vector<Receiver>& receivers() const;
    std::vector<Receiver>& receivers();
    int32 num_receivers() const;
    const Receiver& receiver(int32 idx) const;
    Receiver& receiver(int32 idx);

    void add_shot(const Shot& shot);
    void add_receiver(const Receiver& receiver);

    void clear_shots();
    void clear_receivers();
    void clear();

    bool load_from_segy(const std::string& segy_file);

    void create_simple_survey(
        float64 source_x_start, float64 source_x_end, int32 num_sources,
        float64 receiver_x_start, float64 receiver_x_end, int32 num_receivers,
        float64 depth = 0.0, float64 y = 0.0
    );

private:
    std::vector<Shot> shots_;
    std::vector<Receiver> receivers_;

    float64 apply_scalar(int32 value, int16 scalar) const;
};

}
