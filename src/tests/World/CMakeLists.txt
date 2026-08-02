add_executable(patch_field_diffusion_test
    patch_field_diffusion_test.cpp
)

target_include_directories(patch_field_diffusion_test PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Chemistry
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Diffusion3D/core
    ${CMAKE_CURRENT_SOURCE_DIR}/../../Diffusion3D/core/gpu
)

target_compile_definitions(patch_field_diffusion_test PRIVATE
    IVDBM_CONFIG_DIR="${CMAKE_SOURCE_DIR}/configFiles")

target_link_libraries(patch_field_diffusion_test PRIVATE
    Catch2::Catch2
    Chemistry
    diffusion3d_core_lib
)

set_target_properties(patch_field_diffusion_test PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

add_test(NAME patch_field_diffusion COMMAND patch_field_diffusion_test)

set_tests_properties(patch_field_diffusion
    PROPERTIES
        LABELS "world;cpu;fast"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
