#include "fwi/segy/segy_parser.h"
#include "fwi/geometry/acquisition.h"
#include "fwi/common/ibm_float.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <iomanip>

void create_sample_segy(const std::string& filename) {
    std::cout << "Creating sample SEGY file: " << filename << std::endl;

    std::ofstream file(filename, std::ios::binary);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot create file: " + filename);
    }

    char text_header[3200];
    std::memset(text_header, ' ', 3200);

    std::vector<std::string> text_lines = {
        "C 1 SEISMIC FWI ENGINE - SAMPLE SEGY FILE",
        "C 2 THIS IS A SYNTHETIC FILE FOR DEMONSTRATION",
        "C 3 CREATED FOR TESTING SEGY PARSING FUNCTIONALITY",
        "C 4 DATA FORMAT: IBM FLOAT (FORMAT CODE 1)",
        "C 5 SAMPLE INTERVAL: 1000 MICROSECONDS (1 MS)",
        "C 6 SAMPLES PER TRACE: 250",
        "C 7 ACQUISITION: MARINE STREAMER",
        "C 8 SOURCE TYPE: AIR GUN ARRAY",
        "C 9 WATER DEPTH: 100 METERS",
        "C10 PROCESSING: NONE - RAW FIELD DATA",
        "C11 GEOMETRY: 5 SHOTS, 20 RECEIVERS PER SHOT",
        "C12 COORDINATE SYSTEM: UTM, METERS",
        "C13 SCALAR FOR COORDINATES: -10 (DIVIDE BY 10)",
        "C14 END OF TEXTUAL HEADER"
    };

    for (size_t i = 0; i < text_lines.size() && i < 40; i++) {
        std::string line = text_lines[i];
        if (line.size() > 80) line = line.substr(0, 80);
        while (line.size() < 80) line += ' ';
        std::memcpy(text_header + i * 80, line.c_str(), 80);
    }

    file.write(text_header, 3200);

    char binary_header[400];
    std::memset(binary_header, 0, 400);

    int32_t job_id = 20240001;
    int32_t line_number = 100;
    int32_t reel_number = 1;
    int16_t sample_interval = 1000;
    int16_t samples_per_trace = 250;
    int16_t data_format = 1;
    int16_t ensemble_fold = 20;
    int16_t measurement_system = 1;

    std::memcpy(binary_header + 0, &job_id, 4);
    std::memcpy(binary_header + 4, &line_number, 4);
    std::memcpy(binary_header + 8, &reel_number, 4);
    std::memcpy(binary_header + 16, &sample_interval, 2);
    std::memcpy(binary_header + 20, &samples_per_trace, 2);
    std::memcpy(binary_header + 24, &data_format, 2);
    std::memcpy(binary_header + 26, &ensemble_fold, 2);
    std::memcpy(binary_header + 54, &measurement_system, 2);

    file.write(binary_header, 400);

    int num_shots = 5;
    int num_receivers_per_shot = 20;
    int trace_num = 0;

    for (int shot = 0; shot < num_shots; shot++) {
        double source_x = 1000.0 + shot * 200.0;
        double source_z = 10.0;

        for (int rec = 0; rec < num_receivers_per_shot; rec++) {
            trace_num++;

            double receiver_x = source_x + 50.0 + rec * 25.0;
            double receiver_z = 5.0;

            char trace_header[240];
            std::memset(trace_header, 0, 240);

            int32_t seq_line = trace_num;
            int32_t seq_file = trace_num;
            int32_t field_record = shot + 1;
            int32_t trace_num_field = rec + 1;
            int32_t cdp_num = trace_num;
            int16_t trace_id = 1;
            int32_t sx = static_cast<int32_t>(source_x * 10);
            int32_t gx = static_cast<int32_t>(receiver_x * 10);
            int16_t scalar_coord = -10;
            int16_t scalar_elev = -10;
            int32_t depth_source = static_cast<int32_t>(source_z * 10);
            int32_t rec_elev = static_cast<int32_t>(receiver_z * 10);
            int16_t tr_samples = 250;
            int16_t tr_dt = 1000;

            std::memcpy(trace_header + 0, &seq_line, 4);
            std::memcpy(trace_header + 4, &seq_file, 4);
            std::memcpy(trace_header + 8, &field_record, 4);
            std::memcpy(trace_header + 12, &trace_num_field, 4);
            std::memcpy(trace_header + 20, &cdp_num, 4);
            std::memcpy(trace_header + 28, &trace_id, 2);
            std::memcpy(trace_header + 72, &sx, 4);
            std::memcpy(trace_header + 76, &sx, 4);
            std::memcpy(trace_header + 80, &gx, 4);
            std::memcpy(trace_header + 84, &gx, 4);
            std::memcpy(trace_header + 88, &depth_source, 4);
            std::memcpy(trace_header + 92, &rec_elev, 4);
            std::memcpy(trace_header + 100, &depth_source, 4);
            std::memcpy(trace_header + 104, &rec_elev, 4);
            std::memcpy(trace_header + 108, &scalar_elev, 2);
            std::memcpy(trace_header + 110, &scalar_coord, 2);
            std::memcpy(trace_header + 114, &tr_samples, 2);
            std::memcpy(trace_header + 116, &tr_dt, 2);

            file.write(trace_header, 240);

            double dt = 0.001;
            double offset = std::abs(receiver_x - source_x);
            double velocity = 1500.0;
            double t0 = offset / velocity;

            for (int i = 0; i < 250; i++) {
                double t = i * dt;
                double value = 0.0;

                if (t > t0 - 0.05 && t < t0 + 0.15) {
                    double tau = (t - t0 - 0.05) / 0.05;
                    value = std::exp(-tau * tau) * std::sin(2.0 * M_PI * 30.0 * (t - t0));
                }

                value += (std::rand() / double(RAND_MAX) - 0.5) * 0.05;

                uint32_t ibm_value = fwi::ieee_to_ibm(static_cast<float>(value));
                file.write(reinterpret_cast<char*>(&ibm_value), 4);
            }
        }
    }

    file.close();
    std::cout << "  Created " << trace_num << " traces" << std::endl;
    std::cout << "  File size: " << 3200 + 400 + trace_num * (240 + 250 * 4) << " bytes" << std::endl;
}

int main() {
    std::cout << "=========================================" << std::endl;
    std::cout << "  SEGY File Reading Example             " << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << std::endl;

    std::string filename = "sample_data.sgy";

    try {
        create_sample_segy(filename);
        std::cout << std::endl;
    } catch (const std::exception& e) {
        std::cout << "ERROR creating sample file: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "Parsing SEGY file..." << std::endl;

    fwi::SegyParser parser;
    if (!parser.open(filename)) {
        std::cout << "ERROR: Cannot open SEGY file" << std::endl;
        return 1;
    }

    if (!parser.parse_headers()) {
        std::cout << "ERROR: Cannot parse headers" << std::endl;
        return 1;
    }

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Text Header (first 10 lines):" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    for (int i = 1; i <= 10; i++) {
        std::string line = parser.text_header().get_line(i);
        std::cout << "  " << std::setw(2) << i << ": " << line << std::endl;
    }

    std::cout << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Binary Header:" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    const auto& bh = parser.binary_header();
    std::cout << "  Job ID:                     " << bh.job_id << std::endl;
    std::cout << "  Line number:                " << bh.line_number << std::endl;
    std::cout << "  Sample interval:            " << bh.sample_interval << " us" << std::endl;
    std::cout << "  Samples per trace:          " << bh.samples_per_trace << std::endl;
    std::cout << "  Data sample format:         " << bh.data_sample_format << " (IBM Float)" << std::endl;
    std::cout << "  Ensemble fold:              " << bh.ensemble_fold << std::endl;
    std::cout << "  Measurement system:         " << (bh.measurement_system == 1 ? "Meters" : "Feet") << std::endl;

    std::cout << "  Data format enum:           " << static_cast<int>(parser.data_format()) << std::endl;

    std::cout << std::endl;
    std::cout << "Parsing all traces..." << std::endl;
    if (!parser.parse_all_traces()) {
        std::cout << "ERROR: Cannot parse traces" << std::endl;
        return 1;
    }

    std::cout << "  Parsed " << parser.num_traces() << " traces" << std::endl;
    std::cout << std::endl;

    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Trace Headers Summary:" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

    const auto& traces = parser.traces();

    std::cout << std::setw(6) << "Trace"
              << std::setw(8) << "Shot"
              << std::setw(12) << "Source X"
              << std::setw(12) << "Group X"
              << std::setw(10) << "Offset"
              << std::setw(10) << "Samples" << std::endl;
    std::cout << std::string(60, '-') << std::endl;

    for (int i = 0; i < std::min(10, (int)traces.size()); i++) {
        const auto& th = traces[i].header;
        double sx = th.source_x / 10.0;
        double gx = th.group_x / 10.0;
        double offset = std::abs(gx - sx);

        std::cout << std::setw(6) << i + 1
                  << std::setw(8) << th.field_record_number
                  << std::setw(12) << std::fixed << std::setprecision(1) << sx
                  << std::setw(12) << std::fixed << std::setprecision(1) << gx
                  << std::setw(10) << std::fixed << std::setprecision(1) << offset
                  << std::setw(10) << traces[i].data.size() << std::endl;
    }

    if (traces.size() > 10) {
        std::cout << "  ... and " << traces.size() - 10 << " more traces" << std::endl;
    }

    std::cout << std::endl;
    std::cout << "-----------------------------------------" << std::endl;
    std::cout << "Trace Data Sample (Trace 1):" << std::endl;
    std::cout << "-----------------------------------------" << std::endl;

    if (!traces.empty() && !traces[0].data.empty()) {
        std::cout << "  Sample count: " << traces[0].data.size() << std::endl;
        std::cout << "  First 20 samples:" << std::endl;

        double max_amp = 0.0;
        int max_idx = 0;
        for (int i = 0; i < (int)traces[0].data.size(); i++) {
            if (std::abs(traces[0].data[i]) > max_amp) {
                max_amp = std::abs(traces[0].data[i]);
                max_idx = i;
            }
        }

        std::cout << "  Max amplitude: " << max_amp << " at sample " << max_idx
                  << " (t=" << max_idx * 0.001 << "s)" << std::endl;
        std::cout << std::endl;

        for (int i = 0; i < 20; i++) {
            if (i > 0 && i % 5 == 0) std::cout << std::endl;
            std::cout << "  " << std::fixed << std::setprecision(6) << traces[0].data[i];
        }
        std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Loading geometry from SEGY..." << std::endl;
    fwi::AcquisitionGeometry geometry;
    if (geometry.load_from_segy(filename)) {
        std::cout << "  Loaded " << geometry.num_shots() << " shots and "
                  << geometry.num_receivers() << " receivers" << std::endl;

        if (geometry.num_shots() > 0) {
            std::cout << "  First shot: (" << geometry.shot(0).x()
                      << ", " << geometry.shot(0).z() << ")" << std::endl;
        }
        if (geometry.num_receivers() > 0) {
            std::cout << "  First receiver: (" << geometry.receiver(0).x()
                      << ", " << geometry.receiver(0).z() << ")" << std::endl;
        }
    }

    parser.close();

    std::cout << std::endl;
    std::cout << "=========================================" << std::endl;
    std::cout << "  SEGY parsing completed successfully!   " << std::endl;
    std::cout << "=========================================" << std::endl;

    return 0;
}
