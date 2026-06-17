add_executable(example_forward_modeling
    example_forward_modeling.cpp
)
target_link_libraries(example_forward_modeling PRIVATE
    fwi_core fwi_segy fwi_geometry fwi_modeling fwi_io fwi_common
)
target_compile_definitions(example_forward_modeling PRIVATE _USE_MATH_DEFINES NOMINMAX)

add_executable(example_segy_read
    example_segy_read.cpp
)
target_link_libraries(example_segy_read PRIVATE
    fwi_segy fwi_geometry fwi_common
)
target_compile_definitions(example_segy_read PRIVATE _USE_MATH_DEFINES NOMINMAX)

add_executable(example_velocity_model
    example_velocity_model.cpp
)
target_link_libraries(example_velocity_model PRIVATE
    fwi_modeling fwi_io fwi_common
)
target_compile_definitions(example_velocity_model PRIVATE _USE_MATH_DEFINES NOMINMAX)
