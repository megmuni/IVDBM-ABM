#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "species_registry.h"

#include <cmath>

TEST_CASE("SpeciesRegistry ivdbm_default diffusivity scales with Q", "[chemistry][registry]")
{
    const double Q = 0.85;
    SpeciesRegistry registry = SpeciesRegistry::ivdbm_default(Q);

    REQUIRE(registry.diffusing_species().size() == 3);

    const double d_tnf = registry.diffusivity(0);
    const double d_tgf = registry.diffusivity(1);
    const double d_il1 = registry.diffusivity(2);

    CHECK(d_tnf == Approx(0.0018 * 0.1 * Q));
    CHECK(d_tgf == Approx(0.00156 * 0.1 * Q));
    CHECK(d_il1 == Approx(0.0018 * 0.1 * Q));

    const SpeciesDescriptor &tnf = registry.descriptor(0);
    CHECK(tnf.concentration_channel == 0);
    CHECK(tnf.diffused_channel == 3);
}

TEST_CASE("SpeciesRegistry rejects invalid swelling ratio", "[chemistry][registry]")
{
    SpeciesRegistry registry;
    REQUIRE_THROWS_AS(registry.set_swelling_ratio(0.0), std::invalid_argument);
}
