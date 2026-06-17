#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

#include "fwi/common/types.h"
#include "fwi/segy/segy_parser.h"
#include "fwi/geometry/shot.h"
#include "fwi/geometry/receiver.h"
#include "fwi/geometry/acquisition.h"
#include "fwi/modeling/velocity_model.h"
#include "fwi/modeling/wavelet.h"
#include "fwi/modeling/source.h"
#include "fwi/modeling/acoustic_2d_fd.h"
#include "fwi/io/vtk_writer.h"
#include "fwi/fwi_engine.h"

namespace py = pybind11;

PYBIND11_MODULE(pyfwi, m) {
    m.doc() = "Seismic Full Waveform Inversion (FWI) Engine Python Bindings";

    py::enum_<fwi::SegyDataFormat>(m, "SegyDataFormat")
        .value("IBM_FLOAT", fwi::SegyDataFormat::IBM_FLOAT)
        .value("INT32", fwi::SegyDataFormat::INT32)
        .value("INT16", fwi::SegyDataFormat::INT16)
        .value("IEEE_FLOAT", fwi::SegyDataFormat::IEEE_FLOAT)
        .value("IEEE_DOUBLE", fwi::SegyDataFormat::IEEE_DOUBLE)
        .value("INT24", fwi::SegyDataFormat::INT24)
        .value("INT8", fwi::SegyDataFormat::INT8);

    py::enum_<fwi::WaveletType>(m, "WaveletType")
        .value("RICKER", fwi::WaveletType::RICKER)
        .value("GAUSSIAN", fwi::WaveletType::GAUSSIAN)
        .value("ORMSBY", fwi::WaveletType::ORMSBY)
        .value("SINE", fwi::WaveletType::SINE);

    py::enum_<fwi::SourceType>(m, "SourceType")
        .value("PRESSURE", fwi::SourceType::PRESSURE)
        .value("PARTICLE_VELOCITY_X", fwi::SourceType::PARTICLE_VELOCITY_X)
        .value("PARTICLE_VELOCITY_Z", fwi::SourceType::PARTICLE_VELOCITY_Z)
        .value("EXPLOSIVE", fwi::SourceType::EXPLOSIVE);

    py::enum_<fwi::VtkFormat>(m, "VtkFormat")
        .value("LEGACY_ASCII", fwi::VtkFormat::LEGACY_ASCII)
        .value("LEGACY_BINARY", fwi::VtkFormat::LEGACY_BINARY)
        .value("XML_IMAGE_DATA", fwi::VtkFormat::XML_IMAGE_DATA);

    py::enum_<fwi::ComputeDevice>(m, "ComputeDevice")
        .value("CPU", fwi::ComputeDevice::CPU)
        .value("CUDA", fwi::ComputeDevice::CUDA);

    py::class_<fwi::Point2D>(m, "Point2D")
        .def(py::init<>())
        .def(py::init<double, double>())
        .def_readwrite("x", &fwi::Point2D::x)
        .def_readwrite("z", &fwi::Point2D::z);

    py::class_<fwi::Point3D>(m, "Point3D")
        .def(py::init<>())
        .def(py::init<double, double, double>())
        .def_readwrite("x", &fwi::Point3D::x)
        .def_readwrite("y", &fwi::Point3D::y)
        .def_readwrite("z", &fwi::Point3D::z);

    py::class_<fwi::Grid2D>(m, "Grid2D")
        .def(py::init<>())
        .def(py::init<int, int, double, double, double, double>(),
             py::arg("nx"), py::arg("nz"), py::arg("dx"), py::arg("dz"),
             py::arg("ox") = 0.0, py::arg("oz") = 0.0)
        .def_readwrite("nx", &fwi::Grid2D::nx)
        .def_readwrite("nz", &fwi::Grid2D::nz)
        .def_readwrite("dx", &fwi::Grid2D::dx)
        .def_readwrite("dz", &fwi::Grid2D::dz)
        .def_readwrite("ox", &fwi::Grid2D::ox)
        .def_readwrite("oz", &fwi::Grid2D::oz)
        .def("index", &fwi::Grid2D::index)
        .def("coord_x", &fwi::Grid2D::coord_x)
        .def("coord_z", &fwi::Grid2D::coord_z);

    py::class_<fwi::Shot>(m, "Shot")
        .def(py::init<>())
        .def(py::init<int, const fwi::Point3D&>())
        .def(py::init<int, double, double, double>())
        .def("id", &fwi::Shot::id)
        .def("location", &fwi::Shot::location, py::return_value_policy::reference)
        .def("x", &fwi::Shot::x)
        .def("y", &fwi::Shot::y)
        .def("z", &fwi::Shot::z)
        .def("set_id", &fwi::Shot::set_id)
        .def("set_location", py::overload_cast<const fwi::Point3D&>(&fwi::Shot::set_location))
        .def("set_location", py::overload_cast<double, double, double>(&fwi::Shot::set_location));

    py::class_<fwi::Receiver>(m, "Receiver")
        .def(py::init<>())
        .def(py::init<int, const fwi::Point3D&>())
        .def(py::init<int, double, double, double>())
        .def("id", &fwi::Receiver::id)
        .def("location", &fwi::Receiver::location, py::return_value_policy::reference)
        .def("x", &fwi::Receiver::x)
        .def("y", &fwi::Receiver::y)
        .def("z", &fwi::Receiver::z)
        .def("set_id", &fwi::Receiver::set_id)
        .def("set_location", py::overload_cast<const fwi::Point3D&>(&fwi::Receiver::set_location))
        .def("set_location", py::overload_cast<double, double, double>(&fwi::Receiver::set_location))
        .def("data", py::overload_cast<>(&fwi::Receiver::data), py::return_value_policy::reference)
        .def("set_data", &fwi::Receiver::set_data);

    py::class_<fwi::Wavelet>(m, "Wavelet")
        .def(py::init<>())
        .def_static("ricker", &fwi::Wavelet::ricker,
                    py::arg("dt"), py::arg("peak_freq"),
                    py::arg("t0") = 0.0, py::arg("n_cycles") = 6)
        .def_static("gaussian", &fwi::Wavelet::gaussian,
                    py::arg("dt"), py::arg("std_dev"),
                    py::arg("t0") = 0.0, py::arg("duration") = 2.0)
        .def_static("ormsby", &fwi::Wavelet::ormsby,
                    py::arg("dt"), py::arg("f1"), py::arg("f2"), py::arg("f3"), py::arg("f4"),
                    py::arg("t0") = 0.0, py::arg("duration") = 2.0)
        .def_static("sine", &fwi::Wavelet::sine,
                    py::arg("dt"), py::arg("freq"),
                    py::arg("t0") = 0.0, py::arg("duration") = 2.0)
        .def("type", &fwi::Wavelet::type)
        .def("dt", &fwi::Wavelet::dt)
        .def("num_samples", &fwi::Wavelet::num_samples)
        .def("duration", &fwi::Wavelet::duration)
        .def("peak_frequency", &fwi::Wavelet::peak_frequency)
        .def("t0", &fwi::Wavelet::t0)
        .def("data", py::overload_cast<>(&fwi::Wavelet::data), py::return_value_policy::reference)
        .def("at", &fwi::Wavelet::at)
        .def("at_time", &fwi::Wavelet::at_time)
        .def("normalize", &fwi::Wavelet::normalize, py::arg("max_amp") = 1.0);

    py::class_<fwi::Source>(m, "Source")
        .def(py::init<>())
        .def(py::init<const fwi::Shot&, const fwi::Wavelet&, fwi::SourceType>(),
             py::arg("shot"), py::arg("wavelet"), py::arg("type") = fwi::SourceType::EXPLOSIVE)
        .def("shot", py::overload_cast<>(&fwi::Source::shot), py::return_value_policy::reference)
        .def("wavelet", py::overload_cast<>(&fwi::Source::wavelet), py::return_value_policy::reference)
        .def("type", &fwi::Source::type)
        .def("set_type", &fwi::Source::set_type)
        .def("location", &fwi::Source::location)
        .def("x", &fwi::Source::x)
        .def("y", &fwi::Source::y)
        .def("z", &fwi::Source::z)
        .def("at", &fwi::Source::at)
        .def("at_time", &fwi::Source::at_time);

    py::class_<fwi::VelocityModel2D>(m, "VelocityModel2D")
        .def(py::init<>())
        .def(py::init<const fwi::Grid2D&>())
        .def(py::init<int, int, double, double, double, double>(),
             py::arg("nx"), py::arg("nz"), py::arg("dx"), py::arg("dz"),
             py::arg("ox") = 0.0, py::arg("oz") = 0.0)
        .def("grid", py::overload_cast<>(&fwi::VelocityModel2D::grid), py::return_value_policy::reference)
        .def("velocity", py::overload_cast<>(&fwi::VelocityModel2D::velocity), py::return_value_policy::reference)
        .def("density", py::overload_cast<>(&fwi::VelocityModel2D::density), py::return_value_policy::reference)
        .def("get_velocity", &fwi::VelocityModel2D::get_velocity)
        .def("get_density", &fwi::VelocityModel2D::get_density)
        .def("set_velocity", &fwi::VelocityModel2D::set_velocity)
        .def("set_density", &fwi::VelocityModel2D::set_density)
        .def("set_uniform_velocity", &fwi::VelocityModel2D::set_uniform_velocity)
        .def("set_uniform_density", &fwi::VelocityModel2D::set_uniform_density)
        .def("set_gradient_velocity", &fwi::VelocityModel2D::set_gradient_velocity)
        .def("set_layered_velocity", &fwi::VelocityModel2D::set_layered_velocity)
        .def("add_rectangular_anomaly", &fwi::VelocityModel2D::add_rectangular_anomaly)
        .def("get_max_velocity", &fwi::VelocityModel2D::get_max_velocity)
        .def("get_min_velocity", &fwi::VelocityModel2D::get_min_velocity)
        .def("compute_stable_dt", &fwi::VelocityModel2D::compute_stable_dt, py::arg("cfl") = 0.5)
        .def("compute_pml_size", &fwi::VelocityModel2D::compute_pml_size, py::arg("freq_max"), py::arg("npml") = 20);

    py::class_<fwi::FDSimulationParams>(m, "FDSimulationParams")
        .def(py::init<>())
        .def_readwrite("dt", &fwi::FDSimulationParams::dt)
        .def_readwrite("num_time_steps", &fwi::FDSimulationParams::num_time_steps)
        .def_readwrite("pml_width", &fwi::FDSimulationParams::pml_width)
        .def_readwrite("snapshot_interval", &fwi::FDSimulationParams::snapshot_interval)
        .def_readwrite("record_interval", &fwi::FDSimulationParams::record_interval)
        .def_readwrite("use_pml", &fwi::FDSimulationParams::use_pml)
        .def_readwrite("record_receivers", &fwi::FDSimulationParams::record_receivers);

    py::class_<fwi::Snapshot2D>(m, "Snapshot2D")
        .def(py::init<>())
        .def_readwrite("time_step", &fwi::Snapshot2D::time_step)
        .def_readwrite("time", &fwi::Snapshot2D::time)
        .def_readwrite("pressure", &fwi::Snapshot2D::pressure)
        .def_readwrite("velocity_x", &fwi::Snapshot2D::velocity_x)
        .def_readwrite("velocity_z", &fwi::Snapshot2D::velocity_z);

    py::class_<fwi::AcquisitionGeometry>(m, "AcquisitionGeometry")
        .def(py::init<>())
        .def("shots", py::overload_cast<>(&fwi::AcquisitionGeometry::shots), py::return_value_policy::reference)
        .def("receivers", py::overload_cast<>(&fwi::AcquisitionGeometry::receivers), py::return_value_policy::reference)
        .def("num_shots", &fwi::AcquisitionGeometry::num_shots)
        .def("num_receivers", &fwi::AcquisitionGeometry::num_receivers)
        .def("shot", py::overload_cast<int>(&fwi::AcquisitionGeometry::shot), py::return_value_policy::reference)
        .def("receiver", py::overload_cast<int>(&fwi::AcquisitionGeometry::receiver), py::return_value_policy::reference)
        .def("add_shot", &fwi::AcquisitionGeometry::add_shot)
        .def("add_receiver", &fwi::AcquisitionGeometry::add_receiver)
        .def("clear_shots", &fwi::AcquisitionGeometry::clear_shots)
        .def("clear_receivers", &fwi::AcquisitionGeometry::clear_receivers)
        .def("clear", &fwi::AcquisitionGeometry::clear)
        .def("load_from_segy", &fwi::AcquisitionGeometry::load_from_segy)
        .def("create_simple_survey", &fwi::AcquisitionGeometry::create_simple_survey,
             py::arg("source_x_start"), py::arg("source_x_end"), py::arg("num_sources"),
             py::arg("receiver_x_start"), py::arg("receiver_x_end"), py::arg("num_receivers"),
             py::arg("depth") = 0.0, py::arg("y") = 0.0);

    py::class_<fwi::VtkWriter>(m, "VtkWriter")
        .def(py::init<>())
        .def("set_format", &fwi::VtkWriter::set_format)
        .def("format", &fwi::VtkWriter::format)
        .def("write_snapshot", &fwi::VtkWriter::write_snapshot)
        .def("write_velocity_model", &fwi::VtkWriter::write_velocity_model)
        .def("write_snapshot_series", &fwi::VtkWriter::write_snapshot_series)
        .def("write_3d_snapshot", &fwi::VtkWriter::write_3d_snapshot,
             py::arg("filename"), py::arg("snapshot"), py::arg("grid"),
             py::arg("y_thickness") = 10.0);

    py::class_<fwi::Acoustic2DFD>(m, "Acoustic2DFD")
        .def(py::init<>())
        .def(py::init<const fwi::VelocityModel2D&, const fwi::FDSimulationParams&>())
        .def("set_model", &fwi::Acoustic2DFD::set_model)
        .def("set_params", &fwi::Acoustic2DFD::set_params)
        .def("model", &fwi::Acoustic2DFD::model, py::return_value_policy::reference)
        .def("params", &fwi::Acoustic2DFD::params, py::return_value_policy::reference)
        .def("snapshots", &fwi::Acoustic2DFD::snapshots, py::return_value_policy::reference)
        .def("receiver_data", &fwi::Acoustic2DFD::receiver_data, py::return_value_policy::reference)
        .def("time_axis", &fwi::Acoustic2DFD::time_axis, py::return_value_policy::reference)
        .def("add_source", &fwi::Acoustic2DFD::add_source)
        .def("set_receivers", &fwi::Acoustic2DFD::set_receivers)
        .def("initialize", &fwi::Acoustic2DFD::initialize)
        .def("run", &fwi::Acoustic2DFD::run)
        .def("step", &fwi::Acoustic2DFD::step, py::arg("num_steps") = 1)
        .def("current_pressure", &fwi::Acoustic2DFD::current_pressure, py::return_value_policy::reference)
        .def("current_velocity_x", &fwi::Acoustic2DFD::current_velocity_x, py::return_value_policy::reference)
        .def("current_velocity_z", &fwi::Acoustic2DFD::current_velocity_z, py::return_value_policy::reference)
        .def("current_time", &fwi::Acoustic2DFD::current_time)
        .def("current_time_step", &fwi::Acoustic2DFD::current_time_step);

    py::class_<fwi::FwiEngine>(m, "FwiEngine")
        .def(py::init<>())
        .def("set_compute_device", &fwi::FwiEngine::set_compute_device)
        .def("compute_device", &fwi::FwiEngine::compute_device)
        .def("load_segy_data", &fwi::FwiEngine::load_segy_data)
        .def("create_synthetic_velocity_model", &fwi::FwiEngine::create_synthetic_velocity_model,
             py::arg("nx"), py::arg("nz"), py::arg("dx"), py::arg("dz"),
             py::arg("water_velocity") = 1500.0,
             py::arg("sediment_velocity") = 2500.0,
             py::arg("basement_velocity") = 3500.0)
        .def("set_velocity_model", &fwi::FwiEngine::set_velocity_model)
        .def("velocity_model", &fwi::FwiEngine::velocity_model, py::return_value_policy::reference)
        .def("create_survey", &fwi::FwiEngine::create_survey,
             py::arg("source_x_start"), py::arg("source_x_end"), py::arg("num_sources"),
             py::arg("receiver_x_start"), py::arg("receiver_x_end"), py::arg("num_receivers"),
             py::arg("depth") = 0.0)
        .def("geometry", py::overload_cast<>(&fwi::FwiEngine::geometry), py::return_value_policy::reference)
        .def("set_simulation_params", &fwi::FwiEngine::set_simulation_params)
        .def("simulation_params", &fwi::FwiEngine::simulation_params, py::return_value_policy::reference)
        .def("set_wavelet", &fwi::FwiEngine::set_wavelet)
        .def("wavelet", &fwi::FwiEngine::wavelet, py::return_value_policy::reference)
        .def("run_forward_modeling", &fwi::FwiEngine::run_forward_modeling, py::arg("shot_index") = 0)
        .def("snapshots", &fwi::FwiEngine::snapshots, py::return_value_policy::reference)
        .def("receiver_data", &fwi::FwiEngine::receiver_data, py::return_value_policy::reference)
        .def("time_axis", &fwi::FwiEngine::time_axis, py::return_value_policy::reference)
        .def("export_velocity_model", &fwi::FwiEngine::export_velocity_model)
        .def("export_snapshot", &fwi::FwiEngine::export_snapshot)
        .def("export_all_snapshots", &fwi::FwiEngine::export_all_snapshots,
             py::arg("output_dir"), py::arg("prefix") = "snapshot")
        .def("export_receiver_data", &fwi::FwiEngine::export_receiver_data)
        .def("set_vtk_format", &fwi::FwiEngine::set_vtk_format)
        .def("solver", py::overload_cast<>(&fwi::FwiEngine::solver), py::return_value_policy::reference);

    m.def("ibm_to_ieee", &fwi::ibm_to_ieee);
    m.def("ieee_to_ibm", &fwi::ieee_to_ibm);
    m.def("ebcdic_to_ascii", &fwi::ebcdic_to_ascii);

    m.attr("SEGY_TEXT_HEADER_SIZE") = fwi::SEGY_TEXT_HEADER_SIZE;
    m.attr("SEGY_BINARY_HEADER_SIZE") = fwi::SEGY_BINARY_HEADER_SIZE;
    m.attr("SEGY_TRACE_HEADER_SIZE") = fwi::SEGY_TRACE_HEADER_SIZE;

    py::class_<fwi::SegyParser>(m, "SegyParser")
        .def(py::init<>())
        .def("open", &fwi::SegyParser::open)
        .def("close", &fwi::SegyParser::close)
        .def("is_open", &fwi::SegyParser::is_open)
        .def("parse_headers", &fwi::SegyParser::parse_headers)
        .def("parse_all_traces", &fwi::SegyParser::parse_all_traces)
        .def("parse_trace", &fwi::SegyParser::parse_trace)
        .def("text_header", &fwi::SegyParser::text_header, py::return_value_policy::reference)
        .def("binary_header", &fwi::SegyParser::binary_header, py::return_value_policy::reference)
        .def("traces", &fwi::SegyParser::traces, py::return_value_policy::reference)
        .def("num_traces", &fwi::SegyParser::num_traces)
        .def("data_format", &fwi::SegyParser::data_format);
}
