add_executable(chemistry_tests
    species_registry_test.cpp
    chemical_environment_test.cpp
)

target_link_libraries(chemistry_tests PRIVATE Catch2::Catch2 Chemistry FieldVariable)

add_test(NAME chemistry_species_registry COMMAND chemistry_tests "[chemistry][registry]")
add_test(NAME chemistry_environment COMMAND chemistry_tests "[chemistry][environment]")
