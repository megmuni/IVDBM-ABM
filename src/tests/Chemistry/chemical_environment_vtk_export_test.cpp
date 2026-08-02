#include <catch2/catch.hpp>

#include <fstream>
#include <iterator>
#include <string>

#include "chemical_environment.h"
#include "chemical_environment_vtk_export.h"

#include "../../FieldVariable/Usr_FieldVariables/Chemical.h"
#include "../../common.h"
#include "../../enums.h"

#ifndef IVDBM_CHEM_CONFIG_DIR
#define IVDBM_CHEM_CONFIG_DIR "configFiles"
#endif

namespace
{

std::string test_chem_config_path()
{
    return std::string(IVDBM_CHEM_CONFIG_DIR) + "/simulation_config.template.json";
}

void load_test_env(ChemicalEnvironment &env)
{
    env.load_from_config(test_chem_config_path(), 1.0, 1.0);
    env.allocate_channels_from_config();
}

} // namespace

TEST_CASE("ChemicalEnvironment VTK export writes a .vti file", "[chemistry][vtk]")
{
    ChemicalEnvironment env(2, 2, 2, 0.01);
    load_test_env(env);
    env.set_concentration(0, TNF, 1.f);

    const std::string path = "output/ivdbm_chem_vtk_export_test.vti";
    ChemicalEnvironmentVtkExportOptions options;
    options.export_concentrations = true;

    REQUIRE(export_chemical_environment_to_vti(env, path, options));

    std::ifstream file(path);
    REQUIRE(file.good());
    std::string contents((std::istreambuf_iterator<char>(file)),
                         std::istreambuf_iterator<char>());
    REQUIRE(contents.find("VTKFile") != std::string::npos);
}

TEST_CASE("format_chemical_environment_vti_path builds padded filenames",
          "[chemistry][vtk]")
{
    REQUIRE(format_chemical_environment_vti_path("output", "chem", 3) ==
            "output/chem_t00003.vti");
}
