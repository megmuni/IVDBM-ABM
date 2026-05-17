# Core diffusion3d sources and headers (header-only interface library)

set(DIFFUSION3D_CORE_HEADERS
    scalar_field_grid.h
    multi_species_field_grid.h
    multi_species_diffusion_settings.h
    multi_species_diffusion_engine.h
    explicit_multi_species_heat_stepper.h
    make_multi_species_diffusion_engine.h
)

add_library(diffusion3d_core INTERFACE)
target_include_directories(diffusion3d_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

# Also add as OBJECT lib for compilation of .cpp files that use headers
add_library(diffusion3d_core_objs OBJECT
    explicit_multi_species_heat_stepper.cpp
)
target_include_directories(diffusion3d_core_objs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(diffusion3d_core_objs PUBLIC diffusion3d_core)
