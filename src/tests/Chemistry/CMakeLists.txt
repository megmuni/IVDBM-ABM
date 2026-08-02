add_executable(chemistry_tests
    chemical_environment_test.cpp
)

target_compile_definitions(chemistry_tests PRIVATE
    IVDBM_CHEM_CONFIG_DIR="${CMAKE_SOURCE_DIR}/configFiles")

target_link_libraries(chemistry_tests PRIVATE Catch2::Catch2 Chemistry FieldVariable)

add_test(NAME chemistry_environment COMMAND chemistry_tests "[chemistry][environment]")
add_test(NAME chemistry_config_load COMMAND chemistry_tests "[chemistry][environment][config]")

set_tests_properties(
    chemistry_environment
    chemistry_config_load
    PROPERTIES
        LABELS "chemistry;cpu;fast"
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
)
