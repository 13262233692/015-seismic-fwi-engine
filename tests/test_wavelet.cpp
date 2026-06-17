#include "fwi/modeling/wavelet.h"
#include <iostream>
#include <cmath>
#include <vector>

int main() {
    std::cout << "=== Testing Wavelet Generation ===" << std::endl;

    int passed = 0;
    int total = 0;

    std::cout << "\n1. Testing Ricker wavelet..." << std::endl;
    total++;
    try {
        double dt = 0.001;
        double peak_freq = 30.0;
        fwi::Wavelet ricker = fwi::Wavelet::ricker(dt, peak_freq);

        std::cout << "  Ricker wavelet created:" << std::endl;
        std::cout << "    Sample interval: " << ricker.dt() << " s" << std::endl;
        std::cout << "    Number of samples: " << ricker.num_samples() << std::endl;
        std::cout << "    Duration: " << ricker.duration() << " s" << std::endl;
        std::cout << "    Peak frequency: " << ricker.peak_frequency() << " Hz" << std::endl;
        std::cout << "    Type: " << static_cast<int>(ricker.type()) << std::endl;

        if (ricker.num_samples() > 0 && std::abs(ricker.dt() - dt) < 1e-10) {
            std::cout << "  PASS: Ricker wavelet parameters correct" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Ricker wavelet parameters incorrect" << std::endl;
        }

        const std::vector<double>& data = ricker.data();
        double max_val = 0.0;
        int max_idx = 0;
        for (int i = 0; i < data.size(); i++) {
            if (std::abs(data[i]) > max_val) {
                max_val = std::abs(data[i]);
                max_idx = i;
            }
        }
        std::cout << "    Maximum amplitude: " << max_val << " at sample " << max_idx << std::endl;

        ricker.normalize(1.0);
        double max_after_norm = 0.0;
        for (double v : ricker.data()) {
            if (std::abs(v) > max_after_norm) max_after_norm = std::abs(v);
        }
        total++;
        if (std::abs(max_after_norm - 1.0) < 1e-6) {
            std::cout << "  PASS: Normalization works correctly" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Normalization failed, max = " << max_after_norm << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n2. Testing Gaussian wavelet..." << std::endl;
    total++;
    try {
        double dt = 0.001;
        double std_dev = 0.02;
        fwi::Wavelet gaussian = fwi::Wavelet::gaussian(dt, std_dev);

        std::cout << "  Gaussian wavelet created:" << std::endl;
        std::cout << "    Sample interval: " << gaussian.dt() << " s" << std::endl;
        std::cout << "    Number of samples: " << gaussian.num_samples() << std::endl;
        std::cout << "    Duration: " << gaussian.duration() << " s" << std::endl;

        if (gaussian.num_samples() > 0) {
            std::cout << "  PASS: Gaussian wavelet created successfully" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Gaussian wavelet has no samples" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n3. Testing Ormsby wavelet..." << std::endl;
    total++;
    try {
        double dt = 0.001;
        fwi::Wavelet ormsby = fwi::Wavelet::ormsby(dt, 5.0, 10.0, 40.0, 50.0);

        std::cout << "  Ormsby wavelet created:" << std::endl;
        std::cout << "    Sample interval: " << ormsby.dt() << " s" << std::endl;
        std::cout << "    Number of samples: " << ormsby.num_samples() << std::endl;

        if (ormsby.num_samples() > 0) {
            std::cout << "  PASS: Ormsby wavelet created successfully" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Ormsby wavelet has no samples" << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
    }

    std::cout << "\n4. Testing Sine wavelet..." << std::endl;
    total++;
    try {
        double dt = 0.001;
        double freq = 25.0;
        fwi::Wavelet sine = fwi::Wavelet::sine(dt, freq);

        std::cout << "  Sine wavelet created:" << std::endl;
        std::cout << "    Sample interval: " << sine.dt() << " s" << std::endl;
        std::cout << "    Number of samples: " << sine.num_samples() << std::endl;
        std::cout << "    Frequency: " << freq << " Hz" << std::endl;

        if (sine.num_samples() > 0) {
            std::cout << "  PASS: Sine wavelet created successfully" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: Sine wavelet has no samples" << std::endl;
        }

        total++;
        double t = 0.05;
        double expected = std::sin(2.0 * M_PI * freq * t);
        double actual = sine.at_time(t);
        if (std::abs(expected - actual) < 0.1) {
            std::cout << "  PASS: at_time interpolation works" << std::endl;
            passed++;
        } else {
            std::cout << "  FAIL: at_time interpolation failed, expected=" << expected << ", actual=" << actual << std::endl;
        }
    } catch (const std::exception& e) {
        std::cout << "  FAIL: Exception: " << e.what() << std::endl;
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
