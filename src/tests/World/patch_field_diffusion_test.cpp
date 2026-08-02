#define CATCH_CONFIG_MAIN
#include <catch2/catch.hpp>

#include "patch_field_diffusion.h"

#include "../../Chemistry/chemical_environment_config.h"
#include "../../Chemistry/diffusion_algorithm_name.h"
#include "../../Chemistry/species_registry.h"
#include "../../enums.h"

#include <cmath>
#include <map>
#include <string>
#include <vector>

#ifndef IVDBM_CONFIG_DIR
#define IVDBM_CONFIG_DIR "configFiles"
#endif

TEST_CASE("PatchFieldDiffusion advances one species over a tick",
          "[world][patch_field_diffusion]") {
  const int nx = 4, ny = 4, nz = 4;
  const double h = 0.01;
  const std::size_t n = static_cast<std::size_t>(nx) * ny * nz;

  const std::string config_path = std::string(IVDBM_CONFIG_DIR) +
                                  "/simulation_config.template.json";
  const ChemicalEnvironmentConfig cfg =
      load_chemical_environment_config(config_path);
  SpeciesRegistry registry = SpeciesRegistry::from_config(cfg, 1.0, 1.0);
  PatchFieldDiffusion pfd(nx, ny, nz, h);
  pfd.set_registry(registry);

  std::vector<float> conc_tnf(n, 0.f);
  std::vector<float> conc_tgf(n, 0.f);
  std::vector<float> conc_il1(n, 0.f);
  std::vector<float> conc_o2(n, 0.f);
  std::vector<float> diffused_tnf(n, 0.f);
  std::vector<float> diffused_tgf(n, 0.f);
  std::vector<float> diffused_il1(n, 0.f);
  std::vector<float> diffused_o2(n, 0.f);

  const int cx = nx / 2, cy = ny / 2, cz = nz / 2;
  conc_tnf[static_cast<std::size_t>(cx + nx * (cy + ny * cz))] = 1.f;

  std::map<SpeciesId, SpeciesDiffusionBuffers> buffers;
  buffers[TNF] = {conc_tnf.data(), diffused_tnf.data()};
  buffers[TGF] = {conc_tgf.data(), diffused_tgf.data()};
  buffers[IL1beta] = {conc_il1.data(), diffused_il1.data()};
  buffers[o2] = {conc_o2.data(), diffused_o2.data()};

  pfd.diffuse_all_species(buffers, 30.0);

  float max_diffused = 0.f;
  for (float v : diffused_tnf)
    max_diffused = std::max(max_diffused, std::abs(v));
  REQUIRE(max_diffused > 0.f);
}

TEST_CASE("PatchFieldDiffusion default and custom algorithm selection",
          "[world][patch_field_diffusion]") {
  PatchFieldDiffusion pfd(4, 4, 4, 0.01);
  REQUIRE(pfd.configured_algorithm() ==
          patch_field_default_diffusion_algorithm());

  pfd.set_diffusion_algorithm(DiffusionAlgorithm::ExplicitHeatEquation);
  REQUIRE(pfd.configured_algorithm() ==
          DiffusionAlgorithm::ExplicitHeatEquation);
  REQUIRE(std::string(pfd.configured_algorithm_label()) ==
          std::string(diffusion_algorithm_label(
              DiffusionAlgorithm::ExplicitHeatEquation)));

#ifdef DIFFUSION3D_CUDA
  pfd.set_diffusion_algorithm(DiffusionAlgorithm::GpuStencil);
  REQUIRE(pfd.configured_algorithm() == DiffusionAlgorithm::GpuStencil);
#endif
}
