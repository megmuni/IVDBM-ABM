add_executable(chemistry_tests
    species_registry_test.cpp
    chemical_environment_test.cpp
    chemical_environment_vtk_export_test.cpp
)

target_compile_definitions(chemistry_tests PRIVATE
    IVDBM_CHEM_CONFIG_DIR="${CMAKE_SOURCE_DIR}/configFiles")

target_link_libraries(chemistry_tests PRIVATE Catch2::Catch2 Chemistry FieldVariable)

add_test(NAME chemistry_species_registry COMMAND chemistry_tests "[chemistry][registry]")
add_test(NAME chemistry_environment COMMAND chemistry_tests "[chemistry][environment]")
add_test(NAME chemistry_config_load COMMAND chemistry_tests "[chemistry][environment][config]")
add_test(NAME chemistry_vtk_export COMMAND chemistry_tests "[chemistry][vtk]")

set_tests_properties(
    chemistry_species_registry
    chemistry_environment
    chemistry_config_load
    chemistry_vtk_export
    PROPERTIES
        LABELS "chemistry;cpu;fast"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
