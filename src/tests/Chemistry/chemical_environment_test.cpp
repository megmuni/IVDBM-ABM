#include <catch2/catch.hpp>

#include "../../Chemistry/chemical_environment.h"
#include "../../Chemistry/species_registry.h"
#include "../../FieldVariable/Usr_FieldVariables/Chemical.h"
#include "../../common.h"
#include "../../enums.h"

TEST_CASE("ChemicalEnvironment merge updates concentration from delta channel",
          "[chemistry][environment]") {
  const int nx = 2, ny = 2, nz = 2;
  const int num_channels = 7;

  ChemicalEnvironment env(nx, ny, nz);
  env.load_ivdbm_default(1.0);
  env.allocate_channel_storage(num_channels, pcellgrad);

  env.set_concentration(0, TNF, 10.f);
  env.accumulate_secretion(0, TNF, 2.f);

  Chemical world_chem;
  env.merge_and_reset_secretion();
  env.copy_totals_to(world_chem);

#ifndef CALIBRATION
  REQUIRE(env.concentration_at(0, TNF) == Approx(12.f));
#else
  REQUIRE(env.concentration_at(0, TNF) == Approx(2.2f));
#endif
  REQUIRE(env.channels(TNF).secretion_delta[0] == Approx(0.f));
}

TEST_CASE("ChemicalEnvironment agent API reads and accumulates secretion",
          "[chemistry][environment][agent]") {
  ChemicalEnvironment env(2, 2, 2);
  env.load_ivdbm_default(1.0);
  env.allocate_channel_storage(7, pcellgrad);

  env.set_concentration(3, TGF, 7.f);
  REQUIRE(env.concentration_at(3, TGF) == Approx(7.f));
  REQUIRE(env.concentration_at_channel(3, pTGF) == Approx(7.f));

  env.accumulate_secretion(3, TNF, 1.5f);
  REQUIRE(env.channels(TNF).secretion_delta[3] == Approx(1.5f));

  env.set_concentration(2, TGF, 3.5f);
  env.update_chemotaxis_from_species(TGF);
  REQUIRE(env.chemotaxis_at(2) == Approx(3.5f));
}

TEST_CASE("ChemicalEnvironment baseline helpers set concentration and totals",
          "[chemistry][environment][baseline]") {
  ChemicalEnvironment env(2, 2, 1);
  env.load_ivdbm_default(1.0);
  env.allocate_channel_storage(7, pcellgrad);

  env.set_concentration(1, TGF, 5.f);
  REQUIRE(env.concentration_at(1, TGF) == Approx(5.f));

  env.update_chemotaxis_from_species(TGF);
  REQUIRE(env.chemotaxis_at(1) == Approx(5.f));

  env.recompute_world_totals();
  REQUIRE(env.total_tgf() == Approx(5.f));

  Chemical world_chem;
  env.copy_totals_to(world_chem);
  REQUIRE(world_chem.totalTGF == Approx(5.f));
}
