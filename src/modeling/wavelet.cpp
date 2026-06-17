#include "fwi/modeling/wavelet.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace fwi {

namespace {
inline float64 sinc(float64 x) {
    if (std::abs(x) < 1e-10) {
        return 1.0;
    }
    return std::sin(x) / x;
}
}

Wavelet::Wavelet()
    : type_(WaveletType::RICKER)
    , dt_(0.001)
    , t0_(0.0)
    , peak_freq_(30.0)
    , data_() {}

Wavelet Wavelet::ricker(float64 dt, float64 peak_freq, float64 t0, int32 n_cycles) {
    Wavelet w;
    w.type_ = WaveletType::RICKER;
    w.dt_ = dt;
    w.peak_freq_ = peak_freq;
    w.t0_ = t0;

    float64 period = 1.0 / peak_freq;
    float64 duration = n_cycles * period;
    int32 n = static_cast<int32>(std::floor(duration / dt)) + 1;

    if (t0 <= 0.0) {
        w.t0_ = duration / 2.0;
    }

    w.data_.resize(n);

    for (int32 i = 0; i < n; i++) {
        float64 t = i * dt;
        float64 tau = M_PI * peak_freq * (t - w.t0_);
        float64 tau2 = tau * tau;
        w.data_[i] = (1.0 - 2.0 * tau2) * std::exp(-tau2);
    }

    return w;
}

Wavelet Wavelet::gaussian(float64 dt, float64 std_dev, float64 t0, float64 duration) {
    Wavelet w;
    w.type_ = WaveletType::GAUSSIAN;
    w.dt_ = dt;
    w.t0_ = t0;
    w.peak_freq_ = 1.0 / (2.0 * M_PI * std_dev);

    int32 n = static_cast<int32>(std::floor(duration / dt)) + 1;

    if (t0 <= 0.0) {
        w.t0_ = duration / 2.0;
    }

    w.data_.resize(n);

    float64 two_sigma2 = 2.0 * std_dev * std_dev;
    for (int32 i = 0; i < n; i++) {
        float64 t = i * dt;
        float64 diff = t - w.t0_;
        w.data_[i] = std::exp(-diff * diff / two_sigma2);
    }

    return w;
}

Wavelet Wavelet::ormsby(float64 dt, float64 f1, float64 f2, float64 f3, float64 f4, float64 t0, float64 duration) {
    Wavelet w;
    w.type_ = WaveletType::ORMSBY;
    w.dt_ = dt;
    w.t0_ = t0;
    w.peak_freq_ = (f2 + f3) / 2.0;

    int32 n = static_cast<int32>(std::floor(duration / dt)) + 1;

    if (t0 <= 0.0) {
        w.t0_ = duration / 2.0;
    }

    w.data_.resize(n);

    float64 f42 = f4 - f2;
    float64 f31 = f3 - f1;
    float64 f4_2 = f4 * f4;
    float64 f3_2 = f3 * f3;
    float64 f2_2 = f2 * f2;
    float64 f1_2 = f1 * f1;

    for (int32 i = 0; i < n; i++) {
        float64 t = i * dt - w.t0_;

        float64 s4 = sinc(M_PI * f4 * t);
        float64 s3 = sinc(M_PI * f3 * t);
        float64 s2 = sinc(M_PI * f2 * t);
        float64 s1 = sinc(M_PI * f1 * t);

        float64 term1 = (f4_2 / (f4 - f3)) * s4 * s4;
        float64 term2 = (f3_2 / (f4 - f3)) * s3 * s3;
        float64 term3 = (f2_2 / (f2 - f1)) * s2 * s2;
        float64 term4 = (f1_2 / (f2 - f1)) * s1 * s1;

        w.data_[i] = (term1 - term2) - (term3 - term4);
    }

    return w;
}

Wavelet Wavelet::sine(float64 dt, float64 freq, float64 t0, float64 duration) {
    Wavelet w;
    w.type_ = WaveletType::SINE;
    w.dt_ = dt;
    w.t0_ = t0;
    w.peak_freq_ = freq;

    int32 n = static_cast<int32>(std::floor(duration / dt)) + 1;

    if (t0 <= 0.0) {
        w.t0_ = 0.0;
    }

    w.data_.resize(n);

    float64 omega = 2.0 * M_PI * freq;
    for (int32 i = 0; i < n; i++) {
        float64 t = i * dt - w.t0_;
        w.data_[i] = std::sin(omega * t);
    }

    return w;
}

WaveletType Wavelet::type() const { return type_; }
float64 Wavelet::dt() const { return dt_; }
int32 Wavelet::num_samples() const { return static_cast<int32>(data_.size()); }
float64 Wavelet::duration() const { return data_.size() * dt_; }
float64 Wavelet::peak_frequency() const { return peak_freq_; }
float64 Wavelet::t0() const { return t0_; }

const std::vector<float64>& Wavelet::data() const { return data_; }
std::vector<float64>& Wavelet::data() { return data_; }

float64 Wavelet::at(int32 idx) const {
    if (idx < 0 || idx >= static_cast<int32>(data_.size())) {
        return 0.0;
    }
    return data_[idx];
}

float64 Wavelet::at_time(float64 t) const {
    float64 idx_f = (t - t0_ + duration() / 2.0) / dt_;
    int32 idx = static_cast<int32>(std::floor(idx_f));
    float64 frac = idx_f - idx;

    if (idx < 0 || idx + 1 >= static_cast<int32>(data_.size())) {
        return 0.0;
    }

    return data_[idx] * (1.0 - frac) + data_[idx + 1] * frac;
}

void Wavelet::normalize(float64 max_amp) {
    if (data_.empty()) return;

    float64 current_max = 0.0;
    for (float64 v : data_) {
        current_max = std::max(current_max, std::abs(v));
    }

    if (current_max > 0.0) {
        float64 scale = max_amp / current_max;
        for (float64& v : data_) {
            v *= scale;
        }
    }
}

}
