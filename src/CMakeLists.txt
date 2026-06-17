add_library(fwi_common STATIC
    common/ibm_float.cpp
    common/utils.cpp
)

add_library(fwi_segy STATIC
    segy/segy_parser.cpp
    segy/ebcdic.cpp
)
target_link_libraries(fwi_segy PRIVATE fwi_common)

add_library(fwi_geometry STATIC
    geometry/shot.cpp
    geometry/receiver.cpp
    geometry/acquisition.cpp
)
target_link_libraries(fwi_geometry PRIVATE fwi_common)

add_library(fwi_modeling STATIC
    modeling/velocity_model.cpp
    modeling/source.cpp
    modeling/wavelet.cpp
    modeling/acoustic_2d_fd.cpp
    modeling/fwi_gradient.cpp
)
target_link_libraries(fwi_modeling PRIVATE fwi_common fwi_geometry)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(fwi_modeling PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(fwi_modeling PUBLIC FWI_WITH_OPENMP)
endif()

add_library(fwi_io STATIC
    io/vtk_writer.cpp
)
target_link_libraries(fwi_io PRIVATE fwi_common fwi_modeling)

if(FWI_WITH_CUDA)
    add_library(fwi_cuda STATIC
        cuda/acoustic_2d_fd_cuda.cu
    )
    set_target_properties(fwi_cuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        CUDA_ARCHITECTURES "70;75;80;86;89"
    )
    target_link_libraries(fwi_cuda PRIVATE fwi_common)
endif()

add_library(fwi_core STATIC
    fwi_engine.cpp
)
target_link_libraries(fwi_core PRIVATE
    fwi_common
    fwi_segy
    fwi_geometry
    fwi_modeling
    fwi_io
    $<$<BOOL:${FWI_WITH_CUDA}>:$<TARGET_NAME_IF_EXISTS:fwi_cuda>>
)

target_compile_definitions(fwi_common PUBLIC _USE_MATH_DEFINES NOMINMAX)
target_compile_definitions(fwi_segy PUBLIC _USE_MATH_DEFINES NOMINMAX)
target_compile_definitions(fwi_geometry PUBLIC _USE_MATH_DEFINES NOMINMAX)
target_compile_definitions(fwi_modeling PUBLIC _USE_MATH_DEFINES NOMINMAX)
target_compile_definitions(fwi_io PUBLIC _USE_MATH_DEFINES NOMINMAX)
target_compile_definitions(fwi_core PUBLIC _USE_MATH_DEFINES NOMINMAX)

if(FWI_WITH_CUDA)
    target_compile_definitions(fwi_cuda PUBLIC _USE_MATH_DEFINES NOMINMAX)
endif()

if(FWI_WITH_PYTHON)
    find_package(pybind11 QUIET)
    if(pybind11_FOUND)
        pybind11_add_module(pyfwi
            python/bindings.cpp
        )
        target_link_libraries(pyfwi PRIVATE
            fwi_core
            fwi_common
            fwi_segy
            fwi_geometry
            fwi_modeling
            fwi_io
        )
        target_compile_definitions(pyfwi PRIVATE _USE_MATH_DEFINES)
        set_target_properties(pyfwi PROPERTIES
            LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/python
        )
    else()
        message(WARNING "pybind11 not found, skipping Python bindings")
    endif()
endif()
