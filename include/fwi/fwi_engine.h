#pragma once

#include "fwi/common/types.h"
#include "fwi/segy/segy_parser.h"
#include "fwi/geometry/acquisition.h"
#include "fwi/modeling/velocity_model.h"
#include "fwi/modeling/acoustic_2d_fd.h"
#include "fwi/io/vtk_writer.h"
#include <string>
#include <memory>
#include <vector>

namespace fwi {

enum class ComputeDevice {
    CPU,
    CUDA
};

class FwiEngine {
public:
    FwiEngine();
    ~FwiEngine();

    void set_compute_device(ComputeDevice device);
    ComputeDevice compute_device() const;

    bool load_segy_data(const std::string& segy_file);
    void create_synthetic_velocity_model(int32 nx, int32 nz, float64 dx, float64 dz,
                                         float64 water_velocity = 1500.0,
                                         float64 sediment_velocity = 2500.0,
                                         float64 basement_velocity = 3500.0);

    void set_velocity_model(const VelocityModel2D& model);
    const VelocityModel2D& velocity_model() const;

    void create_survey(float64 source_x_start, float64 source_x_end, int32 num_sources,
                       float64 receiver_x_start, float64 receiver_x_end, int32 num_receivers,
                       float64 depth = 0.0);

    const AcquisitionGeometry& geometry() const;
    AcquisitionGeometry& geometry();

    void set_simulation_params(const FDSimulationParams& params);
    const FDSimulationParams& simulation_params() const;

    void set_wavelet(const Wavelet& wavelet);
    const Wavelet& wavelet() const;

    bool run_forward_modeling(int32 shot_index = 0);

    const std::vector<Snapshot2D>& snapshots() const;
    const std::vector<std::vector<float64>>& receiver_data() const;
    const std::vector<float64>& time_axis() const;

    bool export_velocity_model(const std::string& filename);
    bool export_snapshot(int32 snapshot_index, const std::string& filename);
    bool export_all_snapshots(const std::string& output_dir, const std::string& prefix = "snapshot");
    bool export_receiver_data(const std::string& filename);

    void set_vtk_format(VtkFormat format);

    const Acoustic2DFD& solver() const;
    Acoustic2DFD& solver();

private:
    ComputeDevice compute_device_;
    VelocityModel2D velocity_model_;
    AcquisitionGeometry geometry_;
    FDSimulationParams sim_params_;
    Wavelet wavelet_;
    std::unique_ptr<Acoustic2DFD> solver_;
    VtkWriter vtk_writer_;
    VtkFormat vtk_format_;
};

}
