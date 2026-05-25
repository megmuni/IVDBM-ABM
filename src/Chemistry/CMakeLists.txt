add_library(Chemistry
    species_registry.cpp
    chemical_environment.cpp
)

set_target_properties(Chemistry PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(Chemistry PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
