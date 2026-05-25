#include <catch2/catch.hpp>

#include "../../Chemistry/chemical_environment.h"
#include "../../Chemistry/species_registry.h"
#include "../../common.h"
#include "../../enums.h"
#include "../../FieldVariable/Usr_FieldVariables/Chemical.h"

#include <vector>

TEST_CASE("ChemicalEnvironment merge updates concentration from delta channel",
          "[chemistry][environment]")
{
    const int nx = 2, ny = 2, nz = 2;
    const int n = nx * ny * nz;
    const int num_channels = 7;

    std::vector<std::vector<float>> storage(num_channels, std::vector<float>(n, 0.f));
    std::vector<float *> rows(num_channels);
    for (int c = 0; c < num_channels; ++c)
        rows[c] = storage[c].data();

    ChemicalEnvironment env(nx, ny, nz);
    env.load_ivdbm_default(1.0);
    env.bind_legacy_grid(rows.data(), pcellgrad);

    rows[pTNF][0] = 10.f;
    rows[dTNF][0] = 2.f;

    Chemical world_chem;
    world_chem.pTNF = rows[pTNF];
    world_chem.dTNF = rows[dTNF];
    world_chem.pTGF = rows[pTGF];
    world_chem.dTGF = rows[dTGF];
    world_chem.pIL1beta = rows[pIL1beta];
    world_chem.dIL1beta = rows[dIL1beta];
    world_chem.pcellgrad = rows[pcellgrad];

    env.merge_and_reset_secretion(world_chem);

#ifndef CALIBRATION
    REQUIRE(rows[pTNF][0] == Approx(12.f));
#else
    REQUIRE(rows[pTNF][0] == Approx(2.2f));
#endif
    REQUIRE(rows[dTNF][0] == Approx(0.f));
}

TEST_CASE("ChemicalEnvironment agent API reads and accumulates secretion",
          "[chemistry][environment][agent]")
{
    const int n = 8;
    std::vector<std::vector<float>> storage(7, std::vector<float>(n, 0.f));
    std::vector<float *> rows(7);
    for (int c = 0; c < 7; ++c)
        rows[c] = storage[c].data();

    ChemicalEnvironment env(2, 2, 2);
    env.load_ivdbm_default(1.0);
    env.bind_legacy_grid(rows.data(), pcellgrad);

    rows[pTGF][3] = 7.f;
    REQUIRE(env.concentration_at(3, TGF) == Approx(7.f));
    REQUIRE(env.concentration_at_channel(3, pTGF) == Approx(7.f));

    env.accumulate_secretion(3, TNF, 1.5f);
    REQUIRE(rows[dTNF][3] == Approx(1.5f));
}
