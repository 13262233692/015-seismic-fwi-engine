#pragma once

#include "fwi/common/types.h"
#include <vector>
#include <string>
#include <fstream>
#include <array>

namespace fwi {

constexpr int32 SEGY_TEXT_HEADER_SIZE = 3200;
constexpr int32 SEGY_BINARY_HEADER_SIZE = 400;
constexpr int32 SEGY_TRACE_HEADER_SIZE = 240;
constexpr int32 SEGY_TEXT_HEADER_LINES = 40;
constexpr int32 SEGY_TEXT_HEADER_LINE_LEN = 80;

struct SegyTextHeader {
    std::array<unsigned char, SEGY_TEXT_HEADER_SIZE> raw_data;
    std::string ascii_text;

    void parse();
    std::string get_line(int32 line_num) const;
};

struct SegyBinaryHeader {
    std::array<unsigned char, SEGY_BINARY_HEADER_SIZE> raw_data;

    int32 job_id;
    int32 line_number;
    int32 reel_number;
    int16 data_traces_per_ensemble;
    int16 auxiliary_traces_per_ensemble;
    int16 sample_interval;
    int16 sample_interval_original;
    int16 samples_per_trace;
    int16 samples_per_trace_original;
    int16 data_sample_format;
    int16 ensemble_fold;
    int16 trace_sorting_code;
    int16 vertical_sum_code;
    int16 sweep_frequency_start;
    int16 sweep_frequency_end;
    int16 sweep_length;
    int16 sweep_type_code;
    int16 sweep_trace_number;
    int16 sweep_taper_start;
    int16 sweep_taper_end;
    int16 taper_type;
    int16 correlated_data_flag;
    int16 binary_gain_recovered;
    int16 amplitude_recovery_method;
    int16 measurement_system;
    int16 impulse_signal_polarity;
    int32 time_basis_code;
    int16 format_revision_number;
    int16 fixed_length_trace_flag;
    int16 extended_header_count;

    void parse();
};

struct SegyTraceHeader {
    std::array<unsigned char, SEGY_TRACE_HEADER_SIZE> raw_data;

    int32 trace_sequence_line;
    int32 trace_sequence_file;
    int32 field_record_number;
    int32 trace_number_field_record;
    int32 energy_source_point_number;
    int32 cdp_number;
    int32 cdp_trace_number;
    int16 trace_identification_code;
    int16 trace_summation_count;
    int16 trace_stacking_count;
    int32 inline_number;
    int32 crossline_number;
    int32 shot_point_number;
    int16 shot_point_scalar;
    int16 trace_value_scalar;
    int32 source_x;
    int32 source_y;
    int32 group_x;
    int32 group_y;
    int32 coordinate_depth;
    int32 receiver_group_elevation;
    int32 surface_elevation_source;
    int32 source_depth_below_surface;
    int32 datum_elevation_receiver;
    int32 datum_elevation_source;
    int16 water_depth_source;
    int16 water_depth_group;
    int16 scalar_elevations;
    int16 scalar_coordinates;
    int16 source_measurement_unit;
    int16 group_measurement_unit;
    int16 static_correction_applied;
    int16 total_static_applied;
    int16 lag_time_a;
    int16 lag_time_b;
    int16 delay_recording_time;
    int16 mute_time_start;
    int16 mute_time_end;
    int16 samples_per_trace;
    int16 sample_interval;

    void parse();
};

struct SegyTrace {
    SegyTraceHeader header;
    std::vector<float64> data;
};

enum class SegyDataFormat : int16 {
    IBM_FLOAT = 1,
    INT32 = 2,
    INT16 = 3,
    IEEE_FLOAT = 5,
    IEEE_DOUBLE = 6,
    INT24 = 8,
    INT8 = 9
};

class SegyParser {
public:
    SegyParser();
    ~SegyParser();

    bool open(const std::string& filename);
    void close();

    bool is_open() const;

    bool parse_headers();
    bool parse_all_traces();
    bool parse_trace(int32 trace_index, SegyTrace& trace);

    const SegyTextHeader& text_header() const;
    const SegyBinaryHeader& binary_header() const;
    const std::vector<SegyTrace>& traces() const;
    int32 num_traces() const;

    SegyDataFormat data_format() const;

private:
    std::ifstream file_;
    bool is_open_;
    int64 file_size_;

    SegyTextHeader text_header_;
    SegyBinaryHeader binary_header_;
    std::vector<SegyTrace> traces_;

    bool parse_data_sample(const unsigned char* buf, float64& value, SegyDataFormat format);
};

}
