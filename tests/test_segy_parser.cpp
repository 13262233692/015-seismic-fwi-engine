#include "fwi/segy/segy_parser.h"
#include "fwi/common/ibm_float.h"
#include "fwi/common/utils.h"
#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>

int main() {
    std::cout << "=== Testing SEGY Parser ===" << std::endl;

    int passed = 0;
    int total = 0;

    std::cout << "\n1. Testing creation of synthetic SEGY file..." << std::endl;
    total++;

    std::string test_filename = "test_data.sgy";

    {
        std::ofstream file(test_filename, std::ios::binary);
        if (!file.is_open()) {
            std::cout << "  FAIL: Cannot create test file" << std::endl;
            return 1;
        }

        unsigned char text_header[3200];
        std::memset(text_header, ' ', 3200);
        const char* label = "C 1 TEST SEGY FILE FOR FWI ENGINE";
        std::memcpy(text_header, label, std::strlen(label));
        for (int i = 0; i < 3200; i++) {
            if (text_header[i] == ' ') {
                text_header[i] = 0x40;
            }
        }
        file.write(reinterpret_cast<char*>(text_header), 3200);

        unsigned char binary_header[400];
        std::memset(binary_header, 0, 400);

        int32_t job_id = 12345;
        int16_t sample_interval = 1000;
        int16_t samples_per_trace = 100;
        int16_t data_format = 1;

        fwi::write_big_endian_int32(binary_header + 0, job_id);
        fwi::write_big_endian_int16(binary_header + 16, sample_interval);
        fwi::write_big_endian_int16(binary_header + 20, samples_per_trace);
        fwi::write_big_endian_int16(binary_header + 24, data_format);

        file.write(reinterpret_cast<char*>(binary_header), 400);

        for (int trace_num = 0; trace_num < 5; trace_num++) {
            unsigned char trace_header[240];
            std::memset(trace_header, 0, 240);

            int32_t seq_num = trace_num + 1;
            int32_t source_x = 1000 + trace_num * 100;
            int32_t group_x = 2000 + trace_num * 50;
            int16_t scalar = -10;
            int16_t trace_samples = 100;
            int16_t trace_dt = 1000;

            fwi::write_big_endian_int32(trace_header + 0, seq_num);
            fwi::write_big_endian_int32(trace_header + 72, source_x);
            fwi::write_big_endian_int32(trace_header + 80, group_x);
            fwi::write_big_endian_int16(trace_header + 70, scalar);
            fwi::write_big_endian_int16(trace_header + 114, trace_samples);
            fwi::write_big_endian_int16(trace_header + 116, trace_dt);

            file.write(reinterpret_cast<char*>(trace_header), 240);

            for (int i = 0; i < 100; i++) {
                float value = static_cast<float>(std::sin(0.1 * i + trace_num) * (i + 1) / 100.0);
                uint32_t ibm_value = fwi::ieee_to_ibm(value);
                unsigned char buf[4];
                fwi::write_big_endian_uint32(buf, ibm_value);
                file.write(reinterpret_cast<char*>(buf), 4);
            }
        }

        file.close();
        std::cout << "  PASS: Created test SEGY file" << std::endl;
        passed++;
    }

    std::cout << "\n2. Testing SEGY file parsing..." << std::endl;
    total++;
    try {
        fwi::SegyParser parser;
        if (!parser.open(test_filename)) {
            std::cout << "  FAIL: Cannot open SEGY file" << std::endl;
            return 1;
        }

        if (!parser.is_open()) {
            std::cout << "  FAIL: is_open() returned false" << std::endl;
            return 1;
        }

        if (!parser.parse_headers()) {
            std::cout << "  FAIL: Cannot parse headers" << std::endl;
            return 1;
        }

        std::cout << "  Headers parsed successfully:" << std::endl;
        std::cout << "    Binary header job ID: " << parser.binary_header().job_id << std::endl;
        std::cout << "    Sample interval: " << parser.binary_header().sample_interval << " us" << std::endl;
        std::cout << "    Samples per trace: " << parser.binary_header().samples_per_trace << std::endl;
        std::cout << "    Data format code: " << parser.binary_header().data_sample_format << std::endl;

        fwi::SegyDataFormat format = parser.data_format();
        std::cout << "    Data format enum: " << static_cast<int>(format) << std::endl;

        if (parser.binary_header().job_id == 12345 &&
            parser.binary_header().samples_per_trace == 100) {
            std::cout << "  PASS: Binary header values correct" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Binary header values incorrect" << std::endl;
        }

        total++;
        if (!parser.parse_all_traces()) {
            std::cout << "  FAIL: Cannot parse traces" << std::endl;
            return 1;
        }

        std::cout << "  Parsed " << parser.num_traces() << " traces" << std::endl;

        if (parser.num_traces() == 5) {
            std::cout << "  PASS: Correct number of traces" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Incorrect number of traces, got " << parser.num_traces() << std::endl;
        }

        const auto& traces = parser.traces();
        for (int i = 0; i < static_cast<int>(traces.size()); i++) {
            std::cout << "  Trace " << i << ":" << std::endl;
            std::cout << "    Sequence number: " << traces[i].header.trace_sequence_line << std::endl;
            std::cout << "    Source X: " << traces[i].header.source_x << std::endl;
            std::cout << "    Group X: " << traces[i].header.group_x << std::endl;
            std::cout << "    Data size: " << traces[i].data.size() << " samples" << std::endl;
            if (traces[i].data.size() > 0) {
                std::cout << "    First sample: " << traces[i].data[0] << std::endl;
            }
        }

        total++;
        if (traces[0].header.trace_sequence_line == 1 &&
            traces[2].header.trace_sequence_line == 3 &&
            traces[0].data.size() == 100) {
            std::cout << "  PASS: Trace headers and data correct" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Trace headers or data incorrect" << std::endl;
        }

        parser.close();
        if (!parser.is_open()) {
            std::cout << "  PASS: File closed correctly" << std::endl;
        }

    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\n3. Testing individual trace parsing..." << std::endl;
    total++;
    try {
        fwi::SegyParser parser;
        parser.open(test_filename);
        parser.parse_headers();

        fwi::SegyTrace trace;
        if (parser.parse_trace(2, trace)) {
            std::cout << "  Parsed trace 2:" << std::endl;
            std::cout << "    Sequence number: " << trace.header.trace_sequence_line << std::endl;
            std::cout << "    Data samples: " << trace.data.size() << std::endl;

            if (trace.header.trace_sequence_line == 3 && trace.data.size() == 100) {
                std::cout << "  PASS: Individual trace parsing works" << std::endl;
                passed++;
            } else {
                std::cout << "  FAIL: Individual trace parsing incorrect" << std::endl;
            }
        } else {
            std::cout << "  FAIL: Cannot parse individual trace" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::remove(test_filename.c_str());

    std::cout << "\n=== Test Summary ===" << std::endl;
    std::cout << "Passed: " << passed << "/" << total << std::endl;

    if (passed == total) {
        std::cout << "All tests passed!" << std::endl;
        return 0;
    } else {
        std::cout << "Some tests failed!" << std::endl;
        return 1;
    }
}
