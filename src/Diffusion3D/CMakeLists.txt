cmake_minimum_required(VERSION 3.18)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(DIFFUSION3D_CUDA "Enable CUDA GPU diffusion code path" OFF)
option(DIFFUSION3D_TESTS "Build Diffusion3D Catch2 tests (CTest)" ON)
option(DIFFUSION3D_DEMO "Build ParaView demo executable (diffusion3d_demo)" ON)

# RTX 4050 Laptop (Ada) = sm_89. Override: -DDIFFUSION3D_CUDA_ARCHITECTURES="90" or "native"
set(DIFFUSION3D_CUDA_ARCHITECTURES "89" CACHE STRING
    "CUDA architectures for Diffusion3D (CMake format: 89, 89;90, native, all-major)")

if(DIFFUSION3D_CUDA)
    project(Diffusion3D LANGUAGES CXX CUDA)
    find_package(CUDAToolkit REQUIRED)
    message(STATUS "Diffusion3D CUDA architectures: ${DIFFUSION3D_CUDA_ARCHITECTURES}")
else()
    project(Diffusion3D LANGUAGES CXX)
endif()

add_subdirectory(core)

if(DIFFUSION3D_DEMO)
    add_subdirectory(demo)
endif()

if(DIFFUSION3D_TESTS)
    include(FetchContent)

    FetchContent_Declare(
        Catch2
        URL https://github.com/catchorg/Catch2/archive/refs/tags/v2.13.10.tar.gz
    )
    FetchContent_MakeAvailable(Catch2)

    enable_testing()
    add_subdirectory(tests)
endif()
