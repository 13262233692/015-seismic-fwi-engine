#pragma once

#include "fwi/common/types.h"
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

namespace fwi {

bool is_big_endian();

uint16 byteswap16(uint16 value);
uint32 byteswap32(uint32 value);
uint64 byteswap64(uint64 value);

uint16 read_big_endian_uint16(const unsigned char* buf);
int16 read_big_endian_int16(const unsigned char* buf);
uint32 read_big_endian_uint32(const unsigned char* buf);
int32 read_big_endian_int32(const unsigned char* buf);

void write_big_endian_uint16(unsigned char* buf, uint16 value);
void write_big_endian_int16(unsigned char* buf, int16 value);
void write_big_endian_uint32(unsigned char* buf, uint32 value);
void write_big_endian_int32(unsigned char* buf, int32 value);

template <typename T>
T clamp(T value, T min_val, T max_val) {
    return (std::max)(min_val, (std::min)(max_val, value));
}

template <typename T>
int32 nearest_index(T value, T origin, T spacing) {
    return static_cast<int32>(std::round((value - origin) / spacing));
}

}
