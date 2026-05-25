add_executable(chemistry_tests species_registry_test.cpp)

target_link_libraries(chemistry_tests PRIVATE Catch2::Catch2 Chemistry)

add_test(NAME chemistry_species_registry COMMAND chemistry_tests "[chemistry][registry]")
