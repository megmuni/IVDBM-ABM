# Tests for Diffusion3D

add_executable(test_diffusion3d
    test_multi_species_field_grid.cpp
    test_multi_species_diffusion_engine.cpp
)

target_link_libraries(test_diffusion3d PRIVATE diffusion3d_core_objs)

enable_testing()
add_test(NAME test_diffusion3d COMMAND test_diffusion3d)
