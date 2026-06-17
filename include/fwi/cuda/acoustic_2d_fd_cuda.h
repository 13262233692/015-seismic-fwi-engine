#pragma once

#include "fwi/common/types.h"
#include <vector>

namespace fwi {
namespace cuda {

struct CudaSimulationParams {
    int32 nx;
    int32 nz;
    float64 dx;
    float64 dz;
    float64 dt;
    int32 pml_width;
    float64* d_pressure;
    float64* d_velocity_x;
    float64* d_velocity_z;
    float64* d_density;
    float64* d_modulus;
    float64* d_pml_damping_x;
    float64* d_pml_damping_z;
    bool use_pml;
};

void acoustic_2d_fd_update_vx_cuda(
    const CudaSimulationParams& params,
    cudaStream_t stream = nullptr
);

void acoustic_2d_fd_update_vz_cuda(
    const CudaSimulationParams& params,
    cudaStream_t stream = nullptr
);

void acoustic_2d_fd_update_p_cuda(
    const CudaSimulationParams& params,
    cudaStream_t stream = nullptr
);

void apply_pml_damping_cuda(
    const CudaSimulationParams& params,
    bool apply_to_pressure,
    cudaStream_t stream = nullptr
);

void add_source_cuda(
    const CudaSimulationParams& params,
    int32 source_ix,
    int32 source_iz,
    float64 amplitude,
    cudaStream_t stream = nullptr
);

class Acoustic2DFDCuda {
public:
    Acoustic2DFDCuda();
    ~Acoustic2DFDCuda();

    bool initialize(int32 nx, int32 nz, int32 pml_width = 0);
    void cleanup();

    void set_velocity_model(const Array2D<float64>& velocity, const Array2D<float64>& density);
    void set_pml_coefficients(const Array2D<float64>& damp_x, const Array2D<float64>& damp_z);

    void step(
        float64* pressure,
        float64* velocity_x,
        float64* velocity_z,
        int32 num_steps,
        float64 dx,
        float64 dz,
        float64 dt,
        int32 source_ix = -1,
        int32 source_iz = -1,
        const std::vector<float64>& source_amplitudes = {}
    );

    bool is_initialized() const { return is_initialized_; }

private:
    bool is_initialized_;
    int32 nx_;
    int32 nz_;
    int32 nx_total_;
    int32 nz_total_;
    int32 pml_width_;

    float64* d_pressure_;
    float64* d_velocity_x_;
    float64* d_velocity_z_;
    float64* d_density_;
    float64* d_modulus_;
    float64* d_pml_damping_x_;
    float64* d_pml_damping_z_;

    cudaStream_t stream_;
};

}
}
