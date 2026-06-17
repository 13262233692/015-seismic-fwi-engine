#pragma once

#include <cstdint>

namespace fwi {

float ibm_to_ieee(uint32_t ibm);

uint32_t ieee_to_ibm(float ieee);

double ibm32_to_ieee64(uint32_t ibm);

uint32_t ieee64_to_ibm32(double ieee);

}
