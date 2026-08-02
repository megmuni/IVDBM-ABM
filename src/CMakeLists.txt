cmake_minimum_required (VERSION 2.6)
project (Nyvonna_CPU_ABM)

option(BUILD_DIFFUSION3D "Build Diffusion3D module" ON)
option(BUILD_SRC_TESTS "Build Catch2 tests under src/tests" ON)

if(BUILD_DIFFUSION3D)
    set(DIFFUSION3D_TESTS OFF)
    add_subdirectory(Diffusion3D)
endif()

add_subdirectory(Agent)
add_subdirectory(FieldVariable)
add_subdirectory(Chemistry)
add_subdirectory(World)
add_subdirectory(Patch)
add_subdirectory(ECM)
add_subdirectory(ArrayChain)
add_subdirectory(Utilities)

if(BUILD_SRC_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
