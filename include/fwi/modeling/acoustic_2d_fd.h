#pragma once

#include "fwi/common/types.h"
#include "fwi/modeling/velocity_model.h"
#include "fwi/modeling/source.h"
#include "fwi/modeling/wavelet.h"
#include "fwi/geometry/acquisition.h"
#include <vector>
#include <functional>

namespace fwi {

struct FDSimulationParams {
    float64 dt;
    int32 num_time_steps;
    int32 pml_width;
    int32 snapshot_interval;
    int32 record_interval;
    bool use_pml;
    bool record_receivers;

    FDSimulationParams()
        : dt(0.001)
        , num_time_steps(1000)
        , pml_width(20)
        , snapshot_interval(100)
        , record_interval(1)
        , use_pml(true)
        , record_receivers(true)
    {}
};

struct Snapshot2D {
    int32 time_step;
    float64 time;
    Array2D<float64> pressure;
    Array2D<float64> velocity_x;
    Array2D<float64> velocity_z;
};

class Acoustic2DFD {
public:
    Acoustic2DFD();
    Acoustic2DFD(const VelocityModel2D& model, const FDSimulationParams& params);

    void set_model(const VelocityModel2D& model);
    void set_params(const FDSimulationParams& params);

    const VelocityModel2D& model() const;
    const FDSimulationParams& params() const;

    const std::vector<Snapshot2D>& snapshots() const;
    const std::vector<std::vector<float64>>& receiver_data() const;
    const std::vector<float64>& time_axis() const;

    void add_source(const Source& source);
    void set_receivers(const std::vector<Receiver>& receivers);

    void initialize();
    void run();
    void step(int32 num_steps = 1);

    const Array2D<float64>& current_pressure() const;
    const Array2D<float64>& current_velocity_x() const;
    const Array2D<float64>& current_velocity_z() const;

    float64 current_time() const;
    int32 current_time_step() const;

    using SnapshotCallback = std::function<void(const Snapshot2D&)>;
    void set_snapshot_callback(SnapshotCallback callback);

    using StepCallback = std::function<void(int32 step, float64 time)>;
    void set_step_callback(StepCallback callback);

protected:
    void apply_pml_damping(float64& val_x, float64& val_z, int32 ix, int32 iz);
    void init_pml_coefficients();
    void apply_source(int32 time_idx);
    void record_receivers(int32 time_idx);
    void save_snapshot(int32 time_idx);

    void compute_derivatives_x(const Array2D<float64>& field, Array2D<float64>& deriv, int32 order = 4);
    void compute_derivatives_z(const Array2D<float64>& field, Array2D<float64>& deriv, int32 order = 4);

    void update_velocity_x(int32 nx, int32 nz, Array2D<float64>& vx, const Array2D<float64>& p, const Array2D<float64>& density, float64 dt, float64 dx);
    void update_velocity_z(int32 nx, int32 nz, Array2D<float64>& vz, const Array2D<float64>& p, const Array2D<float64>& density, float64 dt, float64 dz);
    void update_pressure(int32 nx, int32 nz, Array2D<float64>& p, const Array2D<float64>& vx, const Array2D<float64>& vz, const Array2D<float64>& modulus, float64 dt, float64 dx, float64 dz);

private:
    VelocityModel2D model_;
    FDSimulationParams params_;

    int32 nx_;
    int32 nz_;
    int32 nx_total_;
    int32 nz_total_;
    float64 dx_;
    float64 dz_;

    Array2D<float64> pressure_;
    Array2D<float64> pressure_prev_;
    Array2D<float64> velocity_x_;
    Array2D<float64> velocity_z_;
    Array2D<float64> density_;
    Array2D<float64> modulus_;

    Array2D<float64> pml_damping_x_;
    Array2D<float64> pml_damping_z_;
    Array2D<float64> pml_psi_x_;
    Array2D<float64> pml_psi_z_;

    std::vector<Source> sources_;
    std::vector<Receiver> receivers_;
    std::vector<std::pair<int32, int32>> receiver_indices_;

    int32 current_step_;
    float64 current_time_;

    std::vector<Snapshot2D> snapshots_;
    std::vector<std::vector<float64>> receiver_data_;
    std::vector<float64> time_axis_;

    SnapshotCallback snapshot_callback_;
    StepCallback step_callback_;

    int32 fd_order_;
};

}
