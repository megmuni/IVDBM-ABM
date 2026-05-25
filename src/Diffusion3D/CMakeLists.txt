cmake_minimum_required(VERSION 3.18)
project(Diffusion3D LANGUAGES CXX)

# Resolve repo root (src/Diffusion3D -> ../../)
get_filename_component(IVDBM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(LEGACY_DIFFUSION3D_SRC "${IVDBM_ROOT}/diffusion3d/src")

option(DIFFUSION3D_CUDA "Enable CUDA GPU diffusion code path in legacy tests" OFF)
option(DIFFUSION3D_TESTS "Build Diffusion3D Catch2 tests (CTest)" ON)

add_subdirectory(core)

if(DIFFUSION3D_TESTS)
    include(FetchContent)

    # Catch2 v2.x, header-only style include: <catch2/catch.hpp>
    FetchContent_Declare(
        Catch2
        URL https://github.com/catchorg/Catch2/archive/refs/tags/v2.13.10.tar.gz
    )
    FetchContent_MakeAvailable(Catch2)

    enable_testing()
    add_subdirectory(tests)
endif()
