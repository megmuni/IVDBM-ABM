cmake_minimum_required(VERSION 3.18)

include(FetchContent)

option(GPU_DIFFUSE "Enable CUDA GPU diffusion code path" OFF)
option(LOGGING "Enable detailed logging output" OFF)
option(DIFFUSION3D_TESTS "Build diffusion3d tests (CTest)" ON)

if(GPU_DIFFUSE)
    project(diffusion3d LANGUAGES CXX CUDA)
else()
    project(diffusion3d LANGUAGES CXX)
endif()

message(STATUS "Configuring diffusion3d project")

if(GPU_DIFFUSE)
    find_package(CUDAToolkit REQUIRED)
endif()

set(DIFFUSION3D_SOURCES
    src/main.cpp
    src/cpu_diffusion3d.cpp
    src/analytic_diffusion3d.cpp
    src/diffusion3d_timestep.cpp
    src/diffusion3d_solver.cpp
    src/diffusion3d_common.h
    src/diffusion3d_cuda_check.h
    src/gpu_diffusion3d.cuh
    src/logger.h
    src/vtk_export.h
    src/convolutionFFT3D_common.h
    src/analytic_diffusion3d.h
)

if(GPU_DIFFUSE)
    list(APPEND DIFFUSION3D_SOURCES
        src/gpu_diffusion3d.cu
        src/convolutionFFT3D.cu
        src/diffusion3d_fft_scratch.cu
    )
else()
    list(APPEND DIFFUSION3D_SOURCES
        src/gpu_diffusion3d_stub.cpp
        src/diffusion3d_fft_scratch_stub.cpp
    )
endif()

add_executable(diffusion3d ${DIFFUSION3D_SOURCES})
target_include_directories(diffusion3d PRIVATE src)
set_target_properties(diffusion3d PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)

if(GPU_DIFFUSE)
    target_compile_definitions(diffusion3d PRIVATE GPU_DIFFUSE=1)
    target_link_libraries(diffusion3d PRIVATE CUDA::cudart CUDA::cufft)
    set_target_properties(diffusion3d PROPERTIES CUDA_ARCHITECTURES native)
endif()

if(LOGGING)
    message(STATUS "Logging is enabled")
    target_compile_definitions(diffusion3d PRIVATE LOGGING=1)
endif()

if(DIFFUSION3D_TESTS)
    message(STATUS "Catch2 tests moved to src/Diffusion3D/tests — build from repo root with BUILD_DIFFUSION3D=ON")
endif()