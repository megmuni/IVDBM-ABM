add_library(Chemistry
    species_registry.cpp
    chemical_environment_config.cpp
    chemical_environment.cpp
    patch_field_diffusion.cpp
)

set_target_properties(Chemistry PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
)

target_include_directories(Chemistry PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/third_party
)

if(BUILD_DIFFUSION3D)
    target_link_libraries(Chemistry PUBLIC diffusion3d_core_lib)
endif()
