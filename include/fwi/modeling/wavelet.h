#pragma once

#include "fwi/common/types.h"
#include <vector>

namespace fwi {

enum class WaveletType {
    RICKER,
    GAUSSIAN,
    ORMSBY,
    SINE,
    BERlage
};

class Wavelet {
public:
    Wavelet();

    static Wavelet ricker(float64 dt, float64 peak_freq, float64 t0 = 0.0, int32 n_cycles = 6);
    static Wavelet gaussian(float64 dt, float64 std_dev, float64 t0 = 0.0, float64 duration = 2.0);
    static Wavelet ormsby(float64 dt, float64 f1, float64 f2, float64 f3, float64 f4, float64 t0 = 0.0, float64 duration = 2.0);
    static Wavelet sine(float64 dt, float64 freq, float64 t0 = 0.0, float64 duration = 2.0);

    WaveletType type() const;
    float64 dt() const;
    int32 num_samples() const;
    float64 duration() const;
    float64 peak_frequency() const;
    float64 t0() const;

    const std::vector<float64>& data() const;
    std::vector<float64>& data();

    float64 at(int32 idx) const;
    float64 at_time(float64 t) const;

    void normalize(float64 max_amp = 1.0);

private:
    WaveletType type_;
    float64 dt_;
    float64 t0_;
    float64 peak_freq_;
    std::vector<float64> data_;
};

}
