#pragma once

#include "fwi/common/types.h"
#include "fwi/modeling/velocity_model.h"
#include "fwi/modeling/acoustic_2d_fd.h"
#include "fwi/modeling/source.h"
#include "fwi/geometry/acquisition.h"
#include <vector>
#include <functional>

namespace fwi {

struct FwiGradientResult {
    Array2D<float64> velocity_gradient;
    double residual_norm;
    int32 num_time_steps;
};

enum class GradientComputationMode {
    SERIAL,
    OPENMP_UNSAFE,
    OPENMP_ATOMIC,
    OPENMP_REDUCTION
};

class FwiGradient {
public:
    FwiGradient();

    void set_velocity_model(const VelocityModel2D& model);
    void set_params(const FDSimulationParams& params);
    void set_wavelet(const Wavelet& wavelet);
    void set_geometry(const AcquisitionGeometry& geometry);

    void set_num_threads(int32 num_threads);
    int32 num_threads() const;

    void set_computation_mode(GradientComputationMode mode);
    GradientComputationMode computation_mode() const;

    FwiGradientResult compute_gradient(
        int32 shot_index,
        const std::vector<std::vector<float64>>& observed_data,
        const std::vector<Snapshot2D>& forward_wavefield = {}
    );

    FwiGradientResult compute_gradient_multishot(
        const std::vector<int32>& shot_indices,
        const std::vector<std::vector<std::vector<float64>>>& observed_data_all_shots
    );

    const Array2D<float64>& last_gradient() const;
    double last_residual() const;

protected:
    std::vector<Snapshot2D> run_forward(int32 shot_index);
    std::vector<std::vector<float64>> run_adjoint(
        int32 shot_index,
        const std::vector<std::vector<float64>>& residual,
        std::vector<Snapshot2D>& adjoint_wavefield_out
    );

    void compute_gradient_serial(
        const std::vector<Snapshot2D>& forward_wavefield,
        const std::vector<Snapshot2D>& adjoint_wavefield,
        Array2D<float64>& gradient
    );

    void compute_gradient_openmp_unsafe(
        const std::vector<Snapshot2D>& forward_wavefield,
        const std::vector<Snapshot2D>& adjoint_wavefield,
        Array2D<float64>& gradient
    );

    void compute_gradient_openmp_atomic(
        const std::vector<Snapshot2D>& forward_wavefield,
        const std::vector<Snapshot2D>& adjoint_wavefield,
        Array2D<float64>& gradient
    );

    void compute_gradient_openmp_reduction(
        const std::vector<Snapshot2D>& forward_wavefield,
        const std::vector<Snapshot2D>& adjoint_wavefield,
        Array2D<float64>& gradient
    );

    void compute_residual(
        const std::vector<std::vector<float64>>& observed,
        const std::vector<std::vector<float64>>& synthetic,
        std::vector<std::vector<float64>>& residual,
        double& residual_norm
    );

private:
    VelocityModel2D model_;
    FDSimulationParams params_;
    Wavelet wavelet_;
    AcquisitionGeometry geometry_;

    int32 num_threads_;
    GradientComputationMode mode_;

    Array2D<float64> last_gradient_;
    double last_residual_;

    int32 nx_;
    int32 nz_;
    int32 nx_total_;
    int32 nz_total_;
    int32 pml_;
    float64 dx_;
    float64 dz_;
    float64 dt_;
};

}
