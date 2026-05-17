# Core diffusion3d sources

set(DIFFUSION3D_CORE_SOURCES
    scalar_field_grid.h
    multi_species_field_grid.h
    multi_species_diffusion_settings.h
    multi_species_diffusion_engine.h
    explicit_multi_species_heat_stepper.h
    make_multi_species_diffusion_engine.h
)

add_library(diffusion3d_core_objs OBJECT
    ${DIFFUSION3D_CORE_SOURCES}
)

target_include_directories(diffusion3d_core_objs PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
