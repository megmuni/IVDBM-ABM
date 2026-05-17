# Core CMakeLists.txt for Diffusion3D

add_subdirectory(core)
add_subdirectory(adapters)
add_subdirectory(tests)

# Main library target
add_library(diffusion3d_core INTERFACE)
target_include_directories(diffusion3d_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/core)
