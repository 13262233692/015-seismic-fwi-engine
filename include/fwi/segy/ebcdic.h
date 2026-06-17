#pragma once

#include "fwi/common/types.h"
#include <string>

namespace fwi {

char ebcdic_to_ascii(unsigned char ebcdic_char);

void ebcdic_to_ascii_string(const unsigned char* ebcdic_buf, size_t len, std::string& ascii_str);

void ascii_to_ebcdic_string(const std::string& ascii_str, unsigned char* ebcdic_buf, size_t len);

}
