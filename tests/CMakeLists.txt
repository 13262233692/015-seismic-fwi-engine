add_executable(test_segy_parser
    test_segy_parser.cpp
)
target_link_libraries(test_segy_parser PRIVATE fwi_segy fwi_common)
target_compile_definitions(test_segy_parser PRIVATE _USE_MATH_DEFINES NOMINMAX)
add_test(NAME test_segy_parser COMMAND test_segy_parser)

add_executable(test_wavelet
    test_wavelet.cpp
)
target_link_libraries(test_wavelet PRIVATE fwi_modeling fwi_common)
target_compile_definitions(test_wavelet PRIVATE _USE_MATH_DEFINES NOMINMAX)
add_test(NAME test_wavelet COMMAND test_wavelet)

add_executable(test_velocity_model
    test_velocity_model.cpp
)
target_link_libraries(test_velocity_model PRIVATE fwi_modeling fwi_common)
target_compile_definitions(test_velocity_model PRIVATE _USE_MATH_DEFINES NOMINMAX)
add_test(NAME test_velocity_model COMMAND test_velocity_model)

add_executable(test_acoustic_2d_fd
    test_acoustic_2d_fd.cpp
)
target_link_libraries(test_acoustic_2d_fd PRIVATE
    fwi_modeling fwi_geometry fwi_io fwi_common
)
target_compile_definitions(test_acoustic_2d_fd PRIVATE _USE_MATH_DEFINES NOMINMAX)
add_test(NAME test_acoustic_2d_fd COMMAND test_acoustic_2d_fd)

add_executable(test_ibm_float
    test_ibm_float.cpp
)
target_link_libraries(test_ibm_float PRIVATE fwi_common)
target_compile_definitions(test_ibm_float PRIVATE _USE_MATH_DEFINES NOMINMAX)
add_test(NAME test_ibm_float COMMAND test_ibm_float)

add_executable(test_vtk_writer
    test_vtk_writer.cpp
)
target_link_libraries(test_vtk_writer PRIVATE fwi_io fwi_modeling fwi_common)
target_compile_definitions(test_vtk_writer PRIVATE _USE_MATH_DEFINES NOMINMAX)
add_test(NAME test_vtk_writer COMMAND test_vtk_writer)

add_executable(test_fwi_gradient
    test_fwi_gradient.cpp
)
target_link_libraries(test_fwi_gradient PRIVATE fwi_core fwi_modeling fwi_geometry fwi_io fwi_segy fwi_common)
target_compile_definitions(test_fwi_gradient PRIVATE _USE_MATH_DEFINES NOMINMAX)
add_test(NAME test_fwi_gradient COMMAND test_fwi_gradient)
