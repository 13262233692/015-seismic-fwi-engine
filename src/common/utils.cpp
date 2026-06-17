#include "fwi/common/utils.h"

namespace fwi {

bool is_big_endian() {
    union {
        uint32 i;
        unsigned char c[4];
    } u = {0x01020304};
    return u.c[0] == 1;
}

uint16 byteswap16(uint16 value) {
    return static_cast<uint16>((value << 8) | (value >> 8));
}

uint32 byteswap32(uint32 value) {
    return ((value & 0x000000FF) << 24) |
           ((value & 0x0000FF00) << 8) |
           ((value & 0x00FF0000) >> 8) |
           ((value & 0xFF000000) >> 24);
}

uint64 byteswap64(uint64 value) {
    return ((value & 0x00000000000000FFULL) << 56) |
           ((value & 0x000000000000FF00ULL) << 40) |
           ((value & 0x0000000000FF0000ULL) << 24) |
           ((value & 0x00000000FF000000ULL) << 8) |
           ((value & 0x000000FF00000000ULL) >> 8) |
           ((value & 0x0000FF0000000000ULL) >> 24) |
           ((value & 0x00FF000000000000ULL) >> 40) |
           ((value & 0xFF00000000000000ULL) >> 56);
}

uint16 read_big_endian_uint16(const unsigned char* buf) {
    return static_cast<uint16>((buf[0] << 8) | buf[1]);
}

int16 read_big_endian_int16(const unsigned char* buf) {
    return static_cast<int16>(read_big_endian_uint16(buf));
}

uint32 read_big_endian_uint32(const unsigned char* buf) {
    return (static_cast<uint32>(buf[0]) << 24) |
           (static_cast<uint32>(buf[1]) << 16) |
           (static_cast<uint32>(buf[2]) << 8) |
           static_cast<uint32>(buf[3]);
}

int32 read_big_endian_int32(const unsigned char* buf) {
    return static_cast<int32>(read_big_endian_uint32(buf));
}

void write_big_endian_uint16(unsigned char* buf, uint16 value) {
    buf[0] = static_cast<unsigned char>((value >> 8) & 0xFF);
    buf[1] = static_cast<unsigned char>(value & 0xFF);
}

void write_big_endian_int16(unsigned char* buf, int16 value) {
    write_big_endian_uint16(buf, static_cast<uint16>(value));
}

void write_big_endian_uint32(unsigned char* buf, uint32 value) {
    buf[0] = static_cast<unsigned char>((value >> 24) & 0xFF);
    buf[1] = static_cast<unsigned char>((value >> 16) & 0xFF);
    buf[2] = static_cast<unsigned char>((value >> 8) & 0xFF);
    buf[3] = static_cast<unsigned char>(value & 0xFF);
}

void write_big_endian_int32(unsigned char* buf, int32 value) {
    write_big_endian_uint32(buf, static_cast<uint32>(value));
}

}
