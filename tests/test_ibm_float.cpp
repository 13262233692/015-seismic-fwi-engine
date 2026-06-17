#include "fwi/common/ibm_float.h"
#include <iostream>
#include <cmath>
#include <vector>

int main() {
    std::cout << "=== Testing IBM Float Conversion ===" << std::endl;

    int passed = 0;
    int total = 0;

    std::cout << "\n1. Testing ibm_to_ieee and ieee_to_ibm round-trip..." << std::endl;
    std::vector<float> test_values = {
        0.0f, 1.0f, -1.0f, 3.14159f, -2.71828f,
        1000.0f, 0.001f, 123.456f, -789.012f,
        1.0e10f, 1.0e-10f
    };

    for (float val : test_values) {
        total++;
        uint32_t ibm = fwi::ieee_to_ibm(val);
        float converted = fwi::ibm_to_ieee(ibm);
        float rel_error = std::abs(val - converted) / (std::abs(val) + 1e-20f);

        if (rel_error < 1e-3f || (val == 0.0f && converted == 0.0f)) {
            std::cout << "  PASS: " << val << " -> ibm(0x" << std::hex << ibm << std::dec << ") -> " << converted
                      << " (rel_error: " << rel_error << ")" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: " << val << " -> ibm(0x" << std::hex << ibm << std::dec << ") -> " << converted
                      << " (rel_error: " << rel_error << ")" << std::endl;
        }
    }

    std::cout << "\n2. Testing ibm32_to_ieee64 and ieee64_to_ibm32 round-trip..." << std::endl;
    std::vector<double> test_values64 = {
        0.0, 1.0, -1.0, 3.14159265358979, -2.718281828459045,
        1000.0, 0.001, 12345.6789, -9876.54321,
        1.0e20, 1.0e-20
    };

    for (double val : test_values64) {
        total++;
        uint32_t ibm = fwi::ieee64_to_ibm32(val);
        double converted = fwi::ibm32_to_ieee64(ibm);
        double rel_error = std::abs(val - converted) / (std::abs(val) + 1e-20);

        if (rel_error < 1e-3 || (val == 0.0 && converted == 0.0)) {
            std::cout << "  PASS: " << val << " -> ibm(0x" << std::hex << ibm << std::dec << ") -> " << converted
                      << " (rel_error: " << rel_error << ")" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: " << val << " -> ibm(0x" << std::hex << ibm << std::dec << ") -> " << converted
                      << " (rel_error: " << rel_error << ")" << std::endl;
        }
    }

    std::cout << "\n3. Testing special values..." << std::endl;

    total++;
    uint32_t zero_ibm = fwi::ieee_to_ibm(0.0f);
    float zero_back = fwi::ibm_to_ieee(zero_ibm);
    if (zero_back == 0.0f) {
        std::cout << "  PASS: Zero conversion works correctly" << std::endl;
        passed++;
    } else {
        std::cout << "  FAIL: Zero conversion failed, got: " << zero_back << std::endl;
    }

    total++;
    float large_val = 1.0e30f;
    uint32_t large_ibm = fwi::ieee_to_ibm(large_val);
    float large_back = fwi::ibm_to_ieee(large_ibm);
    if (large_back > 0.0f && std::isfinite(large_back)) {
        std::cout << "  PASS: Large value handled: " << large_val << " -> " << large_back << std::endl;
        passed++;
    } else {
        std::cout << "  FAIL: Large value handling failed" << std::endl;
    }

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
