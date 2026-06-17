#include "fwi/cuda/acoustic_2d_fd_cuda.h"
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cmath>

namespace fwi {
namespace cuda {

__constant__ double FD_COEFF_STAGGERED[4] = {
    1.0 / 2.0 + 1.0 / 24.0,
    -1.0 / 24.0,
    0.0,
    0.0
};

__global__ void update_vx_kernel(
    int32 nx, int32 nz,
    const double* __restrict__ p,
    double* __restrict__ vx,
    const double* __restrict__ density,
    double dt, double dx
) {
    int32 ix = blockIdx.x * blockDim.x + threadIdx.x;
    int32 iz = blockIdx.y * blockDim.y + threadIdx.y;

    if (ix >= 2 && ix < nx - 2 && iz >= 2 && iz < nz - 2) {
        int64 idx = iz * nx + ix;

        double coeff1 = FD_COEFF_STAGGERED[0];
        double coeff2 = FD_COEFF_STAGGERED[1];
        double dt_over_dx = dt / dx;

        double dp_dx = coeff1 * (p[idx] - p[idx - 1])
                     + coeff2 * (p[idx + 1] - p[idx - 2]);
        dp_dx *= dt_over_dx;

        double rho = density[idx];
        vx[idx] -= dp_dx / rho;
    }
}

__global__ void update_vz_kernel(
    int32 nx, int32 nz,
    const double* __restrict__ p,
    double* __restrict__ vz,
    const double* __restrict__ density,
    double dt, double dz
) {
    int32 ix = blockIdx.x * blockDim.x + threadIdx.x;
    int32 iz = blockIdx.y * blockDim.y + threadIdx.y;

    if (ix >= 2 && ix < nx - 2 && iz >= 2 && iz < nz - 2) {
        int64 idx = iz * nx + ix;

        double coeff1 = FD_COEFF_STAGGERED[0];
        double coeff2 = FD_COEFF_STAGGERED[1];
        double dt_over_dz = dt / dz;

        double dp_dz = coeff1 * (p[idx] - p[idx - nx])
                     + coeff2 * (p[idx + nx] - p[idx - 2 * nx]);
        dp_dz *= dt_over_dz;

        double rho = density[idx];
        vz[idx] -= dp_dz / rho;
    }
}

__global__ void update_p_kernel(
    int32 nx, int32 nz,
    double* __restrict__ p,
    const double* __restrict__ vx,
    const double* __restrict__ vz,
    const double* __restrict__ modulus,
    double dt, double dx, double dz
) {
    int32 ix = blockIdx.x * blockDim.x + threadIdx.x;
    int32 iz = blockIdx.y * blockDim.y + threadIdx.y;

    if (ix >= 2 && ix < nx - 2 && iz >= 2 && iz < nz - 2) {
        int64 idx = iz * nx + ix;

        double coeff1 = FD_COEFF_STAGGERED[0];
        double coeff2 = FD_COEFF_STAGGERED[1];
        double dt_over_dx = dt / dx;
        double dt_over_dz = dt / dz;

        double dvx_dx = coeff1 * (vx[idx + 1] - vx[idx])
                      + coeff2 * (vx[idx + 2] - vx[idx - 1]);
        dvx_dx *= dt_over_dx;

        double dvz_dz = coeff1 * (vz[idx + nx] - vz[idx])
                      + coeff2 * (vz[idx + 2 * nx] - vz[idx - nx]);
        dvz_dz *= dt_over_dz;

        double K = modulus[idx];
        p[idx] -= K * (dvx_dx + dvz_dz);
    }
}

__global__ void pml_damping_kernel(
    int32 nx, int32 nz,
    double* __restrict__ vx,
    double* __restrict__ vz,
    double* __restrict__ p,
    const double* __restrict__ damp_x,
    const double* __restrict__ damp_z,
    double dt,
    bool apply_to_pressure
) {
    int32 ix = blockIdx.x * blockDim.x + threadIdx.x;
    int32 iz = blockIdx.y * blockDim.y + threadIdx.y;

    if (ix >= 0 && ix < nx && iz >= 0 && iz < nz) {
        int64 idx = iz * nx + ix;

        double dx = damp_x[idx];
        double dz = damp_z[idx];
        double exp_dx = exp(-dx * dt);
        double exp_dz = exp(-dz * dt);

        vx[idx] *= exp_dx;
        vz[idx] *= exp_dz;

        if (apply_to_pressure) {
            p[idx] *= exp_dx * exp_dz;
        }
    }
}

__global__ void add_source_kernel(
    int32 nx, int32 nz,
    double* __restrict__ p,
    const double* __restrict__ modulus,
    const double* __restrict__ density,
    int32 source_ix, int32 source_iz,
    double amplitude, double dt
) {
    if (threadIdx.x == 0 && blockIdx.x == 0) {
        int64 idx = source_iz * nx + source_ix;
        double K = modulus[idx];
        double rho = density[idx];
        double source_term = amplitude * K * dt / rho;
        p[idx] += source_term;
    }
}

void acoustic_2d_fd_update_vx_cuda(
    const CudaSimulationParams& params,
    cudaStream_t stream
) {
    dim3 blockDim(16, 16);
    dim3 gridDim((params.nx + blockDim.x - 1) / blockDim.x,
                 (params.nz + blockDim.y - 1) / blockDim.y);

    update_vx_kernel<<<gridDim, blockDim, 0, stream>>>(
        params.nx, params.nz,
        params.d_pressure,
        params.d_velocity_x,
        params.d_density,
        params.dt, params.dx
    );
}

void acoustic_2d_fd_update_vz_cuda(
    const CudaSimulationParams& params,
    cudaStream_t stream
) {
    dim3 blockDim(16, 16);
    dim3 gridDim((params.nx + blockDim.x - 1) / blockDim.x,
                 (params.nz + blockDim.y - 1) / blockDim.y);

    update_vz_kernel<<<gridDim, blockDim, 0, stream>>>(
        params.nx, params.nz,
        params.d_pressure,
        params.d_velocity_z,
        params.d_density,
        params.dt, params.dz
    );
}

void acoustic_2d_fd_update_p_cuda(
    const CudaSimulationParams& params,
    cudaStream_t stream
) {
    dim3 blockDim(16, 16);
    dim3 gridDim((params.nx + blockDim.x - 1) / blockDim.x,
                 (params.nz + blockDim.y - 1) / blockDim.y);

    update_p_kernel<<<gridDim, blockDim, 0, stream>>>(
        params.nx, params.nz,
        params.d_pressure,
        params.d_velocity_x,
        params.d_velocity_z,
        params.d_modulus,
        params.dt, params.dx, params.dz
    );
}

void apply_pml_damping_cuda(
    const CudaSimulationParams& params,
    bool apply_to_pressure,
    cudaStream_t stream
) {
    dim3 blockDim(16, 16);
    dim3 gridDim((params.nx + blockDim.x - 1) / blockDim.x,
                 (params.nz + blockDim.y - 1) / blockDim.y);

    pml_damping_kernel<<<gridDim, blockDim, 0, stream>>>(
        params.nx, params.nz,
        params.d_velocity_x,
        params.d_velocity_z,
        params.d_pressure,
        params.d_pml_damping_x,
        params.d_pml_damping_z,
        params.dt,
        apply_to_pressure
    );
}

void add_source_cuda(
    const CudaSimulationParams& params,
    int32 source_ix,
    int32 source_iz,
    float64 amplitude,
    cudaStream_t stream
) {
    add_source_kernel<<<1, 1, 0, stream>>>(
        params.nx, params.nz,
        params.d_pressure,
        params.d_modulus,
        params.d_density,
        source_ix, source_iz,
        amplitude, params.dt
    );
}

Acoustic2DFDCuda::Acoustic2DFDCuda()
    : is_initialized_(false)
    , nx_(0), nz_(0), nx_total_(0), nz_total_(0), pml_width_(0)
    , d_pressure_(nullptr)
    , d_velocity_x_(nullptr)
    , d_velocity_z_(nullptr)
    , d_density_(nullptr)
    , d_modulus_(nullptr)
    , d_pml_damping_x_(nullptr)
    , d_pml_damping_z_(nullptr)
    , stream_(nullptr)
{}

Acoustic2DFDCuda::~Acoustic2DFDCuda() {
    cleanup();
}

bool Acoustic2DFDCuda::initialize(int32 nx, int32 nz, int32 pml_width) {
    if (is_initialized_) {
        cleanup();
    }

    nx_ = nx;
    nz_ = nz;
    pml_width_ = pml_width;
    nx_total_ = nx + 2 * pml_width;
    nz_total_ = nz + 2 * pml_width;

    int64 total_size = static_cast<int64>(nx_total_) * nz_total_;
    size_t bytes = static_cast<size_t>(total_size) * sizeof(float64);

    cudaMalloc(&d_pressure_, bytes);
    cudaMalloc(&d_velocity_x_, bytes);
    cudaMalloc(&d_velocity_z_, bytes);
    cudaMalloc(&d_density_, bytes);
    cudaMalloc(&d_modulus_, bytes);
    cudaMalloc(&d_pml_damping_x_, bytes);
    cudaMalloc(&d_pml_damping_z_, bytes);

    cudaStreamCreate(&stream_);

    is_initialized_ = true;
    return true;
}

void Acoustic2DFDCuda::cleanup() {
    if (!is_initialized_) return;

    cudaFree(d_pressure_);
    cudaFree(d_velocity_x_);
    cudaFree(d_velocity_z_);
    cudaFree(d_density_);
    cudaFree(d_modulus_);
    cudaFree(d_pml_damping_x_);
    cudaFree(d_pml_damping_z_);

    if (stream_) {
        cudaStreamDestroy(stream_);
    }

    d_pressure_ = nullptr;
    d_velocity_x_ = nullptr;
    d_velocity_z_ = nullptr;
    d_density_ = nullptr;
    d_modulus_ = nullptr;
    d_pml_damping_x_ = nullptr;
    d_pml_damping_z_ = nullptr;
    stream_ = nullptr;

    is_initialized_ = false;
}

void Acoustic2DFDCuda::set_velocity_model(const Array2D<float64>& velocity, const Array2D<float64>& density) {
    if (!is_initialized_) return;

    Array2D<float64> density_full(nx_total_, nz_total_);
    Array2D<float64> modulus_full(nx_total_, nz_total_);

    for (int32 iz = 0; iz < nz_; iz++) {
        for (int32 ix = 0; ix < nx_; ix++) {
            int32 ix_full = ix + pml_width_;
            int32 iz_full = iz + pml_width_;
            float64 vel = velocity(ix, iz);
            float64 rho = density(ix, iz);
            density_full(ix_full, iz_full) = rho;
            modulus_full(ix_full, iz_full) = rho * vel * vel;
        }
    }

    for (int32 iz = 0; iz < nz_total_; iz++) {
        for (int32 ix = 0; ix < nx_total_; ix++) {
            if (density_full(ix, iz) <= 0.0) {
                int32 ix_model = clamp(ix - pml_width_, 0, nx_ - 1);
                int32 iz_model = clamp(iz - pml_width_, 0, nz_ - 1);
                density_full(ix, iz) = density(ix_model, iz_model);
                float64 vel = velocity(ix_model, iz_model);
                modulus_full(ix, iz) = density_full(ix, iz) * vel * vel;
            }
        }
    }

    int64 total_size = static_cast<int64>(nx_total_) * nz_total_;
    size_t bytes = static_cast<size_t>(total_size) * sizeof(float64);

    cudaMemcpy(d_density_, density_full.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_modulus_, modulus_full.data(), bytes, cudaMemcpyHostToDevice);
}

void Acoustic2DFDCuda::set_pml_coefficients(const Array2D<float64>& damp_x, const Array2D<float64>& damp_z) {
    if (!is_initialized_) return;

    int64 total_size = static_cast<int64>(nx_total_) * nz_total_;
    size_t bytes = static_cast<size_t>(total_size) * sizeof(float64);

    cudaMemcpy(d_pml_damping_x_, damp_x.data(), bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_pml_damping_z_, damp_z.data(), bytes, cudaMemcpyHostToDevice);
}

void Acoustic2DFDCuda::step(
    float64* pressure,
    float64* velocity_x,
    float64* velocity_z,
    int32 num_steps,
    float64 dx,
    float64 dz,
    float64 dt,
    int32 source_ix,
    int32 source_iz,
    const std::vector<float64>& source_amplitudes
) {
    if (!is_initialized_) return;

    int64 total_size = static_cast<int64>(nx_total_) * nz_total_;
    size_t bytes = static_cast<size_t>(total_size) * sizeof(float64);

    cudaMemcpy(d_pressure_, pressure, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_velocity_x_, velocity_x, bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_velocity_z_, velocity_z, bytes, cudaMemcpyHostToDevice);

    CudaSimulationParams params;
    params.nx = nx_total_;
    params.nz = nz_total_;
    params.dx = dx;
    params.dz = dz;
    params.dt = dt;
    params.pml_width = pml_width_;
    params.d_pressure = d_pressure_;
    params.d_velocity_x = d_velocity_x_;
    params.d_velocity_z = d_velocity_z_;
    params.d_density = d_density_;
    params.d_modulus = d_modulus_;
    params.d_pml_damping_x = d_pml_damping_x_;
    params.d_pml_damping_z = d_pml_damping_z_;
    params.use_pml = pml_width_ > 0;

    for (int32 step = 0; step < num_steps; step++) {
        if (source_ix >= 0 && source_iz >= 0 && step < static_cast<int32>(source_amplitudes.size())) {
            add_source_cuda(params, source_ix, source_iz, source_amplitudes[step], stream_);
        }

        acoustic_2d_fd_update_vx_cuda(params, stream_);
        acoustic_2d_fd_update_vz_cuda(params, stream_);

        if (params.use_pml) {
            apply_pml_damping_cuda(params, false, stream_);
        }

        acoustic_2d_fd_update_p_cuda(params, stream_);

        if (params.use_pml) {
            apply_pml_damping_cuda(params, true, stream_);
        }
    }

    cudaMemcpy(pressure, d_pressure_, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(velocity_x, d_velocity_x_, bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(velocity_z, d_velocity_z_, bytes, cudaMemcpyDeviceToHost);

    cudaStreamSynchronize(stream_);
}

}
}
