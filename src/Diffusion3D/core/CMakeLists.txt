# Core diffusion3d sources and headers (header-only interface library)

set(DIFFUSION3D_CORE_HEADERS
    scalar_field_grid.h
    multi_species_field_grid.h
    multi_species_diffusion_settings.h
    multi_species_diffusion_engine.h
    explicit_multi_species_heat_stepper.h
    explicit_gpu_stencil_stepper.h
    explicit_gpu_fft_stepper.h
    make_multi_species_diffusion_engine.h
    diffusion_algorithm.h
    diffusion3d_timestep.h
    diffusion3d_step_euler_cpu.h
    vtk_export.h
)

add_library(diffusion3d_core INTERFACE)
target_include_directories(diffusion3d_core INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/gpu
)

add_subdirectory(gpu)

add_library(diffusion3d_core_objs OBJECT
    explicit_multi_species_heat_stepper.cpp
    explicit_gpu_stencil_stepper.cpp
    explicit_gpu_fft_stepper.cpp
    diffusion3d_timestep.cpp
    diffusion3d_step_euler_cpu.cpp
)
target_include_directories(diffusion3d_core_objs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(diffusion3d_core_objs PUBLIC diffusion3d_core diffusion3d_core_gpu)
set_target_properties(diffusion3d_core_objs PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)

# Linkable archive for standalone demo and external consumers (I.8).
add_library(diffusion3d_core_lib STATIC
    $<TARGET_OBJECTS:diffusion3d_core_objs>
    $<TARGET_OBJECTS:diffusion3d_core_gpu>
)
target_link_libraries(diffusion3d_core_lib PUBLIC diffusion3d_core)
if(DIFFUSION3D_CUDA)
    target_link_libraries(diffusion3d_core_lib PUBLIC CUDA::cudart CUDA::cufft)
endif()
set_target_properties(diffusion3d_core_lib PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
)
