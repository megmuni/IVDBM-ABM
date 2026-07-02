cmake_minimum_required(VERSION 3.18)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(DIFFUSION3D_CUDA "Enable CUDA GPU diffusion code path" OFF)
option(DIFFUSION3D_TESTS "Build Diffusion3D Catch2 tests (CTest)" ON)
option(DIFFUSION3D_DEMO "Build ParaView demo executable (diffusion3d_demo)" ON)

# RTX 4050 Laptop (Ada) = sm_89. Override: -DDIFFUSION3D_CUDA_ARCHITECTURES="90" or "native"
set(DIFFUSION3D_CUDA_ARCHITECTURES "89" CACHE STRING
    "CUDA architectures for Diffusion3D (CMake format: 89, 89;90, native, all-major)")

project(Diffusion3D LANGUAGES CXX)

if(DIFFUSION3D_CUDA)
    find_package(CUDAToolkit QUIET)
    if(NOT CUDAToolkit_FOUND)
        message(FATAL_ERROR
            "DIFFUSION3D_CUDA=ON but the CUDA toolkit was not found (nvcc missing).\n"
            "On DRAC / Alliance Canada login nodes, load CUDA before running cmake, e.g.:\n"
            "  module load StdEnv/2023 gcc/12.3 cuda/12.2\n"
            "  # or: module load cuda\n"
            "For CPU-only builds, omit -DDIFFUSION3D_CUDA=ON (default is OFF).")
    endif()
    set(CMAKE_CUDA_COMPILER "${CUDAToolkit_NVCC_EXECUTABLE}" CACHE FILEPATH "" FORCE)
    enable_language(CUDA)
    message(STATUS "Diffusion3D CUDA compiler: ${CMAKE_CUDA_COMPILER}")
    message(STATUS "Diffusion3D CUDA architectures: ${DIFFUSION3D_CUDA_ARCHITECTURES}")
endif()

add_subdirectory(core)

if(DIFFUSION3D_DEMO)
    add_subdirectory(demo)
endif()

# Standalone configure (-S src/Diffusion3D): tests live under ../tests
if(DIFFUSION3D_TESTS)
    enable_testing()
    add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/../tests ${CMAKE_BINARY_DIR}/src_tests)
endif()
