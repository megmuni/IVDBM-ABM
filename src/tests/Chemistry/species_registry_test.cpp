#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "chemical_environment_config.h"
#include "species_registry.h"

#include <cmath>

#ifndef IVDBM_CHEM_CONFIG_DIR
#define IVDBM_CHEM_CONFIG_DIR "configFiles"
#endif

namespace
{

std::string test_chem_config_path()
{
    return std::string(IVDBM_CHEM_CONFIG_DIR) + "/chemical_environment.template.json";
}

} // namespace

TEST_CASE("SpeciesRegistry from_config diffusivity scales with Q",
          "[chemistry][registry]")
{
    const double Q = 0.85;
    const ChemicalEnvironmentConfig cfg =
        load_chemical_environment_config(test_chem_config_path());
    SpeciesRegistry registry = SpeciesRegistry::from_config(cfg, Q, 10.0);

    REQUIRE(registry.diffusing_species().size() == 4);

    const double d_tnf = registry.diffusivity(0);
    const double d_tgf = registry.diffusivity(1);
    const double d_il1 = registry.diffusivity(2);

    CHECK(d_tnf == Approx(0.00018 * Q));
    CHECK(d_tgf == Approx(0.000156 * Q));
    CHECK(d_il1 == Approx(0.00018 * Q));

    const SpeciesDescriptor &tnf = registry.descriptor(0);
    CHECK(tnf.concentration_channel == 0);
    CHECK(tnf.diffused_channel == 4);
}

TEST_CASE("SpeciesRegistry O2 diffusivity uses logarithmic stiffness model",
          "[chemistry][registry]")
{
    const double E = 5.0;
    const ChemicalEnvironmentConfig cfg =
        load_chemical_environment_config(test_chem_config_path());
    SpeciesRegistry registry = SpeciesRegistry::from_config(cfg, 1.0, E);

    const double d_o2 = registry.diffusivity(3);
    CHECK(d_o2 == Approx(-0.002 * std::log(E) + 0.0218));
}

TEST_CASE("SpeciesRegistry rejects invalid swelling ratio", "[chemistry][registry]")
{
    SpeciesRegistry registry;
    REQUIRE_THROWS_AS(registry.set_swelling_ratio(0.0), std::invalid_argument);
}
