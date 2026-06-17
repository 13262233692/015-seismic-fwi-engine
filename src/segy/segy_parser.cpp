#include "fwi/segy/segy_parser.h"
#include "fwi/common/utils.h"
#include "fwi/common/ibm_float.h"
#include "fwi/segy/ebcdic.h"

#include <stdexcept>
#include <cstring>
#include <iostream>

namespace fwi {

void SegyTextHeader::parse() {
    ebcdic_to_ascii_string(raw_data.data(), raw_data.size(), ascii_text);
}

std::string SegyTextHeader::get_line(int32 line_num) const {
    if (line_num < 0 || line_num >= SEGY_TEXT_HEADER_LINES) {
        throw std::out_of_range("Line number out of range");
    }
    int32 start = line_num * SEGY_TEXT_HEADER_LINE_LEN;
    return ascii_text.substr(start, SEGY_TEXT_HEADER_LINE_LEN);
}

void SegyBinaryHeader::parse() {
    const unsigned char* buf = raw_data.data();

    job_id = read_big_endian_int32(buf + 0);
    line_number = read_big_endian_int32(buf + 4);
    reel_number = read_big_endian_int32(buf + 8);
    data_traces_per_ensemble = read_big_endian_int16(buf + 12);
    auxiliary_traces_per_ensemble = read_big_endian_int16(buf + 14);
    sample_interval = read_big_endian_int16(buf + 16);
    sample_interval_original = read_big_endian_int16(buf + 18);
    samples_per_trace = read_big_endian_int16(buf + 20);
    samples_per_trace_original = read_big_endian_int16(buf + 22);
    data_sample_format = read_big_endian_int16(buf + 24);
    ensemble_fold = read_big_endian_int16(buf + 26);
    trace_sorting_code = read_big_endian_int16(buf + 28);
    vertical_sum_code = read_big_endian_int16(buf + 30);
    sweep_frequency_start = read_big_endian_int16(buf + 32);
    sweep_frequency_end = read_big_endian_int16(buf + 34);
    sweep_length = read_big_endian_int16(buf + 36);
    sweep_type_code = read_big_endian_int16(buf + 38);
    sweep_trace_number = read_big_endian_int16(buf + 40);
    sweep_taper_start = read_big_endian_int16(buf + 42);
    sweep_taper_end = read_big_endian_int16(buf + 44);
    taper_type = read_big_endian_int16(buf + 46);
    correlated_data_flag = read_big_endian_int16(buf + 48);
    binary_gain_recovered = read_big_endian_int16(buf + 50);
    amplitude_recovery_method = read_big_endian_int16(buf + 52);
    measurement_system = read_big_endian_int16(buf + 54);
    impulse_signal_polarity = read_big_endian_int16(buf + 56);
    time_basis_code = read_big_endian_int32(buf + 60);
    format_revision_number = read_big_endian_int16(buf + 300);
    fixed_length_trace_flag = read_big_endian_int16(buf + 302);
    extended_header_count = read_big_endian_int16(buf + 304);
}

void SegyTraceHeader::parse() {
    const unsigned char* buf = raw_data.data();

    trace_sequence_line = read_big_endian_int32(buf + 0);
    trace_sequence_file = read_big_endian_int32(buf + 4);
    field_record_number = read_big_endian_int32(buf + 8);
    trace_number_field_record = read_big_endian_int32(buf + 12);
    energy_source_point_number = read_big_endian_int32(buf + 16);
    cdp_number = read_big_endian_int32(buf + 20);
    cdp_trace_number = read_big_endian_int32(buf + 24);
    trace_identification_code = read_big_endian_int16(buf + 28);
    trace_summation_count = read_big_endian_int16(buf + 30);
    trace_stacking_count = read_big_endian_int16(buf + 32);
    inline_number = read_big_endian_int32(buf + 188);
    crossline_number = read_big_endian_int32(buf + 192);
    shot_point_number = read_big_endian_int32(buf + 168);
    shot_point_scalar = read_big_endian_int16(buf + 172);
    trace_value_scalar = read_big_endian_int16(buf + 174);
    source_x = read_big_endian_int32(buf + 72);
    source_y = read_big_endian_int32(buf + 76);
    group_x = read_big_endian_int32(buf + 80);
    group_y = read_big_endian_int32(buf + 84);
    coordinate_depth = read_big_endian_int32(buf + 40);
    receiver_group_elevation = read_big_endian_int32(buf + 44);
    surface_elevation_source = read_big_endian_int32(buf + 48);
    source_depth_below_surface = read_big_endian_int32(buf + 52);
    datum_elevation_receiver = read_big_endian_int32(buf + 56);
    datum_elevation_source = read_big_endian_int32(buf + 60);
    water_depth_source = read_big_endian_int16(buf + 144);
    water_depth_group = read_big_endian_int16(buf + 146);
    scalar_elevations = read_big_endian_int16(buf + 68);
    scalar_coordinates = read_big_endian_int16(buf + 70);
    source_measurement_unit = read_big_endian_int16(buf + 150);
    group_measurement_unit = read_big_endian_int16(buf + 152);
    static_correction_applied = read_big_endian_int16(buf + 104);
    total_static_applied = read_big_endian_int16(buf + 106);
    lag_time_a = read_big_endian_int16(buf + 108);
    lag_time_b = read_big_endian_int16(buf + 110);
    delay_recording_time = read_big_endian_int16(buf + 112);
    mute_time_start = read_big_endian_int16(buf + 114);
    mute_time_end = read_big_endian_int16(buf + 116);
    samples_per_trace = read_big_endian_int16(buf + 118);
    sample_interval = read_big_endian_int16(buf + 120);
}

SegyParser::SegyParser() : is_open_(false), file_size_(0) {}

SegyParser::~SegyParser() {
    close();
}

bool SegyParser::open(const std::string& filename) {
    close();

    file_.open(filename, std::ios::binary | std::ios::ate);
    if (!file_.is_open()) {
        return false;
    }

    file_size_ = file_.tellg();
    file_.seekg(0, std::ios::beg);
    is_open_ = true;

    return true;
}

void SegyParser::close() {
    if (file_.is_open()) {
        file_.close();
    }
    is_open_ = false;
    file_size_ = 0;
    traces_.clear();
}

bool SegyParser::is_open() const {
    return is_open_;
}

bool SegyParser::parse_headers() {
    if (!is_open_) {
        return false;
    }

    file_.seekg(0, std::ios::beg);

    file_.read(reinterpret_cast<char*>(text_header_.raw_data.data()), SEGY_TEXT_HEADER_SIZE);
    if (file_.gcount() != SEGY_TEXT_HEADER_SIZE) {
        return false;
    }
    text_header_.parse();

    file_.read(reinterpret_cast<char*>(binary_header_.raw_data.data()), SEGY_BINARY_HEADER_SIZE);
    if (file_.gcount() != SEGY_BINARY_HEADER_SIZE) {
        return false;
    }
    binary_header_.parse();

    return true;
}

bool SegyParser::parse_data_sample(const unsigned char* buf, float64& value, SegyDataFormat format) {
    switch (format) {
        case SegyDataFormat::IBM_FLOAT: {
            uint32 ibm = read_big_endian_uint32(buf);
            value = ibm32_to_ieee64(ibm);
            return true;
        }
        case SegyDataFormat::IEEE_FLOAT: {
            uint32 ieee_bits = read_big_endian_uint32(buf);
            float32 f;
            std::memcpy(&f, &ieee_bits, sizeof(f));
            value = static_cast<float64>(f);
            return true;
        }
        case SegyDataFormat::IEEE_DOUBLE: {
            uint64 ieee_bits = 0;
            for (int i = 0; i < 8; i++) {
                ieee_bits = (ieee_bits << 8) | buf[i];
            }
            float64 d;
            std::memcpy(&d, &ieee_bits, sizeof(d));
            value = d;
            return true;
        }
        case SegyDataFormat::INT32: {
            value = static_cast<float64>(read_big_endian_int32(buf));
            return true;
        }
        case SegyDataFormat::INT16: {
            value = static_cast<float64>(read_big_endian_int16(buf));
            return true;
        }
        case SegyDataFormat::INT24: {
            int32 v = (static_cast<int32>(buf[0]) << 24) |
                      (static_cast<int32>(buf[1]) << 16) |
                      (static_cast<int32>(buf[2]) << 8);
            v >>= 8;
            value = static_cast<float64>(v);
            return true;
        }
        case SegyDataFormat::INT8: {
            value = static_cast<float64>(static_cast<int8_t>(buf[0]));
            return true;
        }
        default:
            return false;
    }
}

bool SegyParser::parse_trace(int32 trace_index, SegyTrace& trace) {
    if (!is_open_) {
        return false;
    }

    SegyDataFormat format = data_format();
    int32 samples_per_trace = binary_header_.samples_per_trace;

    int64 data_offset = SEGY_TEXT_HEADER_SIZE + SEGY_BINARY_HEADER_SIZE;
    int32 bytes_per_sample = 0;

    switch (format) {
        case SegyDataFormat::IBM_FLOAT:
        case SegyDataFormat::IEEE_FLOAT:
        case SegyDataFormat::INT32:
            bytes_per_sample = 4;
            break;
        case SegyDataFormat::IEEE_DOUBLE:
            bytes_per_sample = 8;
            break;
        case SegyDataFormat::INT16:
            bytes_per_sample = 2;
            break;
        case SegyDataFormat::INT24:
            bytes_per_sample = 3;
            break;
        case SegyDataFormat::INT8:
            bytes_per_sample = 1;
            break;
        default:
            return false;
    }

    int64 trace_size = SEGY_TRACE_HEADER_SIZE + bytes_per_sample * samples_per_trace;
    int64 trace_offset = data_offset + static_cast<int64>(trace_index) * trace_size;

    if (trace_offset + trace_size > file_size_) {
        return false;
    }

    file_.seekg(trace_offset, std::ios::beg);

    file_.read(reinterpret_cast<char*>(trace.header.raw_data.data()), SEGY_TRACE_HEADER_SIZE);
    if (file_.gcount() != SEGY_TRACE_HEADER_SIZE) {
        return false;
    }
    trace.header.parse();

    std::vector<unsigned char> data_buf(static_cast<size_t>(bytes_per_sample) * samples_per_trace);
    file_.read(reinterpret_cast<char*>(data_buf.data()), data_buf.size());
    if (file_.gcount() != static_cast<std::streamsize>(data_buf.size())) {
        return false;
    }

    trace.data.resize(samples_per_trace);
    for (int32 i = 0; i < samples_per_trace; i++) {
        if (!parse_data_sample(data_buf.data() + i * bytes_per_sample, trace.data[i], format)) {
            return false;
        }
    }

    return true;
}

bool SegyParser::parse_all_traces() {
    if (!is_open_) {
        return false;
    }

    if (!parse_headers()) {
        return false;
    }

    SegyDataFormat format = data_format();
    int32 samples_per_trace = binary_header_.samples_per_trace;

    int64 data_offset = SEGY_TEXT_HEADER_SIZE + SEGY_BINARY_HEADER_SIZE;
    int32 bytes_per_sample = 0;

    switch (format) {
        case SegyDataFormat::IBM_FLOAT:
        case SegyDataFormat::IEEE_FLOAT:
        case SegyDataFormat::INT32:
            bytes_per_sample = 4;
            break;
        case SegyDataFormat::IEEE_DOUBLE:
            bytes_per_sample = 8;
            break;
        case SegyDataFormat::INT16:
            bytes_per_sample = 2;
            break;
        case SegyDataFormat::INT24:
            bytes_per_sample = 3;
            break;
        case SegyDataFormat::INT8:
            bytes_per_sample = 1;
            break;
        default:
            return false;
    }

    int64 trace_size = SEGY_TRACE_HEADER_SIZE + bytes_per_sample * samples_per_trace;
    int64 remaining = file_size_ - data_offset;
    int32 num_traces = static_cast<int32>(remaining / trace_size);

    if (num_traces < 0) {
        return false;
    }

    traces_.resize(num_traces);
    for (int32 i = 0; i < num_traces; i++) {
        if (!parse_trace(i, traces_[i])) {
            traces_.clear();
            return false;
        }
    }

    return true;
}

const SegyTextHeader& SegyParser::text_header() const {
    return text_header_;
}

const SegyBinaryHeader& SegyParser::binary_header() const {
    return binary_header_;
}

const std::vector<SegyTrace>& SegyParser::traces() const {
    return traces_;
}

int32 SegyParser::num_traces() const {
    return static_cast<int32>(traces_.size());
}

SegyDataFormat SegyParser::data_format() const {
    return static_cast<SegyDataFormat>(binary_header_.data_sample_format);
}

}
