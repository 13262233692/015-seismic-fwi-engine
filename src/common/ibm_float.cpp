#include "fwi/common/ibm_float.h"

#include <cstring>
#include <cmath>

namespace fwi {

float ibm_to_ieee(uint32_t ibm) {
    uint32_t sign = ibm & 0x80000000;
    int32_t exponent = static_cast<int32_t>((ibm & 0x7F000000) >> 24);
    uint32_t mantissa = ibm & 0x00FFFFFF;

    if (mantissa == 0) {
        return 0.0f;
    }

    int32_t exp_val = exponent - 64;

    double fraction = 0.0;
    double divisor = 16.0;

    for (int32_t i = 0; i < 6; i++) {
        uint32_t hex_digit = (mantissa >> (20 - i * 4)) & 0x0F;
        fraction += static_cast<double>(hex_digit) / divisor;
        divisor *= 16.0;
    }

    double value = fraction;
    if (exp_val >= 0) {
        for (int32_t i = 0; i < exp_val; i++) {
            value *= 16.0;
        }
    } else {
        for (int32_t i = 0; i < -exp_val; i++) {
            value /= 16.0;
        }
    }

    if (sign) {
        value = -value;
    }

    return static_cast<float>(value);
}

uint32_t ieee_to_ibm(float ieee) {
    if (ieee == 0.0f) {
        return 0;
    }

    uint32_t ieee_bits;
    std::memcpy(&ieee_bits, &ieee, sizeof(ieee_bits));
    uint32_t sign = ieee_bits & 0x80000000;

    double abs_val = std::fabs(static_cast<double>(ieee));

    if (!std::isfinite(abs_val)) {
        return sign | 0x7FFFFFFF;
    }

    int32_t ibm_exponent = 64;

    while (abs_val >= 1.0) {
        abs_val /= 16.0;
        ibm_exponent++;
    }
    while (abs_val < 1.0 / 16.0) {
        abs_val *= 16.0;
        ibm_exponent--;
    }

    if (ibm_exponent > 127) {
        return sign | 0x7FFFFFFF;
    }
    if (ibm_exponent < 0) {
        return 0;
    }

    uint32_t mantissa = 0;
    double remaining = abs_val;

    for (int32_t i = 0; i < 6; i++) {
        remaining *= 16.0;
        uint32_t digit = static_cast<uint32_t>(remaining);
        if (digit > 15) digit = 15;
        mantissa = (mantissa << 4) | digit;
        remaining -= digit;
    }

    mantissa &= 0x00FFFFFF;

    return sign | (static_cast<uint32_t>(ibm_exponent) << 24) | mantissa;
}

double ibm32_to_ieee64(uint32_t ibm) {
    return static_cast<double>(ibm_to_ieee(ibm));
}

uint32_t ieee64_to_ibm32(double ieee) {
    return ieee_to_ibm(static_cast<float>(ieee));
}

}
