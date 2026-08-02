/**
 * @file ex_chem_random_spikes.cpp
 * @brief Place random chemical concentration spikes on a patch grid and
 * optionally run diffusion ticks (Diffusion3D via ChemicalEnvironment).
 *
 * Run from the project root so default config paths resolve, e.g.:
 *   ./build/bin/ex_chem_random_spikes --nx 32 --ny 32 --nz 32 --spikes 25
 * --ticks 5 Defaults use low diffusivity scale (0.01) so spikes spread slowly
 * on a 32^3 grid.
 *
 * ParaView: File → Open → output/chem_spikes/chem_spikes_t*.vti → group as time
 * series.
 */

#include "../src/Chemistry/chemical_environment.h"
#include "../src/Chemistry/chemical_environment_vtk_export.h"
#include "../src/Chemistry/diffusion_algorithm_name.h"
#include "../src/Diffusion3D/core/diffusion_algorithm.h"
#include "../src/enums.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#ifndef IVDBM_CHEM_CONFIG_DIR
#define IVDBM_CHEM_CONFIG_DIR "configFiles"
#endif

namespace {

struct Options {
  int nx = 32;
  int ny = 32;
  int nz = 32;
  double h_mm = 0.01;
  int num_spikes = 20;
  int num_ticks = 0;
  float spike_min = 0.1f;
  float spike_max = 10.f;
  unsigned seed = 0;
  bool random_seed = true;
  int species = -1; // -1 = random species per spike
  std::string config_path =
      std::string(IVDBM_CHEM_CONFIG_DIR) + "/chemical_environment.json";
  std::string output_dir = "output/chem_spikes";
  bool write_paraview = true;
  double viz_multiplier = 1.0;
  /** Scales D_eff = base_diffusivity * scale (demo default << 1 for slower
   * spread). */
  double diffusivity_scale = 0.01;
  /** Diffusion macro-step in minutes; 0 = use tick_interval_minutes from JSON.
   */
  double tick_dt_minutes = 0.0;
  /** -1 = patch_field default; else force ExplicitHeatEquation / GpuStencil /
   * GpuFftPrecomputed. */
  int algorithm_choice = -1;
};

void usage(const char *prog) {
  std::cerr
      << "Usage: " << prog << " [options]\n"
      << "  --nx N              Grid size X (default 32)\n"
      << "  --ny N              Grid size Y (default 32)\n"
      << "  --nz N              Grid size Z (default 32)\n"
      << "  --h MM              Patch spacing mm (default 0.01)\n"
      << "  --spikes N          Number of random spikes (default 20)\n"
      << "  --ticks N           Diffusion ticks after spikes (default 0)\n"
      << "  --spike-min V       Min spike amount (default 0.1, typical "
         "production < 10)\n"
      << "  --spike-max V       Max spike amount (default 10)\n"
      << "  --seed N            RNG seed (default: random_device)\n"
      << "  --species NAME      TNF | TGF | IL1beta (default: random per "
         "spike)\n"
      << "  --chem-config PATH  chemical_environment.json\n"
      << "  --out-dir PATH      VTK output directory (default "
         "output/chem_spikes)\n"
      << "  --no-paraview       Skip VTK (.vti) export\n"
      << "  --viz-multiplier M  Extra ParaView array *= M when M != 1 (default "
         "1)\n"
      << "  --diffusivity-scale S  D_eff = base_D * S (default 0.01, slower "
         "spread)\n"
      << "  --tick-dt MIN       Diffusion step length in minutes (default: "
         "JSON)\n"
      << "  --algorithm NAME    cpu | stencil | fft (default: build default)\n";
}

bool parse_int(const char *s, int &out) {
  char *end = nullptr;
  const long v = std::strtol(s, &end, 10);
  if (end == s || *end != '\0')
    return false;
  out = static_cast<int>(v);
  return true;
}

bool parse_float(const char *s, float &out) {
  char *end = nullptr;
  const float v = std::strtof(s, &end);
  if (end == s || *end != '\0')
    return false;
  out = v;
  return true;
}

bool parse_double(const char *s, double &out) {
  char *end = nullptr;
  const double v = std::strtod(s, &end);
  if (end == s || *end != '\0')
    return false;
  out = v;
  return true;
}

bool parse_args(int argc, char **argv, Options &opt) {
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    auto need = [&](const char *name) -> const char * {
      if (i + 1 >= argc) {
        std::cerr << "Missing value for " << name << "\n";
        return nullptr;
      }
      return argv[++i];
    };

    if (arg == "--help" || arg == "-h") {
      usage(argv[0]);
      return false;
    }
    if (arg == "--nx") {
      const char *v = need("--nx");
      if (!v || !parse_int(v, opt.nx))
        return false;
    } else if (arg == "--ny") {
      const char *v = need("--ny");
      if (!v || !parse_int(v, opt.ny))
        return false;
    } else if (arg == "--nz") {
      const char *v = need("--nz");
      if (!v || !parse_int(v, opt.nz))
        return false;
    } else if (arg == "--h") {
      const char *v = need("--h");
      if (!v || !parse_double(v, opt.h_mm))
        return false;
    } else if (arg == "--spikes") {
      const char *v = need("--spikes");
      if (!v || !parse_int(v, opt.num_spikes))
        return false;
    } else if (arg == "--ticks") {
      const char *v = need("--ticks");
      if (!v || !parse_int(v, opt.num_ticks))
        return false;
    } else if (arg == "--spike-min") {
      const char *v = need("--spike-min");
      if (!v || !parse_float(v, opt.spike_min))
        return false;
    } else if (arg == "--spike-max") {
      const char *v = need("--spike-max");
      if (!v || !parse_float(v, opt.spike_max))
        return false;
    } else if (arg == "--seed") {
      const char *v = need("--seed");
      int s = 0;
      if (!v || !parse_int(v, s))
        return false;
      opt.seed = static_cast<unsigned>(s);
      opt.random_seed = false;
    } else if (arg == "--species") {
      const char *v = need("--species");
      if (!v)
        return false;
      const std::string name = v;
      if (name == "TNF")
        opt.species = TNF;
      else if (name == "TGF")
        opt.species = TGF;
      else if (name == "IL1beta" || name == "IL1")
        opt.species = IL1beta;
      else {
        std::cerr << "Unknown species: " << name << "\n";
        return false;
      }
    } else if (arg == "--chem-config") {
      const char *v = need("--chem-config");
      if (!v)
        return false;
      opt.config_path = v;
    } else if (arg == "--out-dir") {
      const char *v = need("--out-dir");
      if (!v)
        return false;
      opt.output_dir = v;
    } else if (arg == "--no-paraview") {
      opt.write_paraview = false;
    } else if (arg == "--viz-multiplier") {
      const char *v = need("--viz-multiplier");
      if (!v || !parse_double(v, opt.viz_multiplier))
        return false;
    } else if (arg == "--diffusivity-scale") {
      const char *v = need("--diffusivity-scale");
      if (!v || !parse_double(v, opt.diffusivity_scale))
        return false;
    } else if (arg == "--tick-dt") {
      const char *v = need("--tick-dt");
      if (!v || !parse_double(v, opt.tick_dt_minutes))
        return false;
    } else if (arg == "--algorithm") {
      const char *v = need("--algorithm");
      if (!v)
        return false;
      if (strcmp(v, "cpu") == 0 || strcmp(v, "explicit") == 0)
        opt.algorithm_choice =
            static_cast<int>(DiffusionAlgorithm::ExplicitHeatEquation);
      else if (strcmp(v, "stencil") == 0 || strcmp(v, "gpu") == 0)
        opt.algorithm_choice = static_cast<int>(DiffusionAlgorithm::GpuStencil);
      else if (strcmp(v, "fft") == 0)
        opt.algorithm_choice =
            static_cast<int>(DiffusionAlgorithm::GpuFftPrecomputed);
      else {
        std::cerr << "Unknown algorithm: " << v
                  << " (use cpu, stencil, or fft)\n";
        return false;
      }
    } else {
      std::cerr << "Unknown option: " << arg << "\n";
      return false;
    }
  }
  return true;
}

int patch_index(int ix, int iy, int iz, int nx, int ny) {
  return ix + nx * (iy + ny * iz);
}

const char *species_name(SpeciesId id) {
  switch (id) {
  case TNF:
    return "TNF";
  case TGF:
    return "TGF";
  case IL1beta:
    return "IL1beta";
  default:
    return "?";
  }
}

void place_random_spikes(ChemicalEnvironment &env, const Options &opt,
                         std::mt19937 &rng) {
  std::uniform_int_distribution<int> dist_x(0, opt.nx - 1);
  std::uniform_int_distribution<int> dist_y(0, opt.ny - 1);
  std::uniform_int_distribution<int> dist_z(0, opt.nz - 1);
  std::uniform_real_distribution<float> dist_amp(opt.spike_min, opt.spike_max);
  std::uniform_int_distribution<int> dist_species(0, 2);

  std::cout << "Placing " << opt.num_spikes << " spikes (seed="
            << (opt.random_seed ? "auto" : std::to_string(opt.seed)) << ")\n";

  for (int i = 0; i < opt.num_spikes; ++i) {
    const int ix = dist_x(rng);
    const int iy = dist_y(rng);
    const int iz = dist_z(rng);
    const int idx = patch_index(ix, iy, iz, opt.nx, opt.ny);
    const SpeciesId sp = (opt.species >= 0)
                             ? static_cast<SpeciesId>(opt.species)
                             : static_cast<SpeciesId>(dist_species(rng));
    const float amp = dist_amp(rng);
    env.set_concentration(idx, sp, amp);

    std::cout << "  spike " << i << ": " << species_name(sp) << " at (" << ix
              << "," << iy << "," << iz << ") idx=" << idx << " amp=" << amp
              << "\n";
  }
}

struct SpeciesStats {
  float max_c = 0.f;
  float sum_c = 0.f;
  int max_ix = 0;
  int max_iy = 0;
  int max_iz = 0;
};

SpeciesStats compute_stats(const ChemicalEnvironment &env, SpeciesId sp, int nx,
                           int ny, int nz) {
  SpeciesStats s;
  for (int iz = 0; iz < nz; ++iz)
    for (int iy = 0; iy < ny; ++iy)
      for (int ix = 0; ix < nx; ++ix) {
        const int idx = patch_index(ix, iy, iz, nx, ny);
        const float c = env.concentration_at(idx, sp);
        s.sum_c += c;
        if (c > s.max_c) {
          s.max_c = c;
          s.max_ix = ix;
          s.max_iy = iy;
          s.max_iz = iz;
        }
      }
  return s;
}

void print_stats(ChemicalEnvironment &env, int nx, int ny, int nz,
                 const char *label) {
  std::cout << label << ":\n";
  for (SpeciesId sp : {TNF, TGF, IL1beta}) {
    const SpeciesStats st = compute_stats(env, sp, nx, ny, nz);
    std::cout << "  " << species_name(sp) << "  max=" << st.max_c << " at ("
              << st.max_ix << "," << st.max_iy << "," << st.max_iz << ")"
              << "  sum=" << st.sum_c << "\n";
  }
  env.recompute_world_totals();
  std::cout << "  world totals: TNF=" << env.total_tnf()
            << " TGF=" << env.total_tgf() << " IL1beta=" << env.total_il1beta()
            << "\n";
}

void write_paraview_timestep(const ChemicalEnvironment &env, const Options &opt,
                             int step) {
  ChemicalEnvironmentVtkExportOptions vtk_opts;
  vtk_opts.viz_multiplier = opt.viz_multiplier;
  const std::string path =
      format_chemical_environment_vti_path(opt.output_dir, "chem_spikes", step);
  if (!export_chemical_environment_to_vti(env, path, vtk_opts)) {
    std::cerr << "VTK export failed: " << path << "\n";
    return;
  }
  std::cout << "ParaView: wrote " << path << "\n";
}

} // namespace

int main(int argc, char **argv) {
  Options opt;
  if (!parse_args(argc, argv, opt))
    return 1;

  if (opt.nx < 1 || opt.ny < 1 || opt.nz < 1) {
    std::cerr << "Grid dimensions must be positive.\n";
    return 1;
  }
  if (opt.num_spikes < 0 || opt.num_ticks < 0) {
    std::cerr << "spikes and ticks must be non-negative.\n";
    return 1;
  }
  if (opt.spike_max < opt.spike_min) {
    std::cerr << "spike-max must be >= spike-min.\n";
    return 1;
  }
  if (opt.h_mm <= 0.0) {
    std::cerr << "Patch spacing --h must be > 0 for diffusion.\n";
    return 1;
  }
  if (opt.diffusivity_scale <= 0.0) {
    std::cerr << "--diffusivity-scale must be > 0.\n";
    return 1;
  }
  if (opt.tick_dt_minutes < 0.0) {
    std::cerr << "--tick-dt must be >= 0.\n";
    return 1;
  }

  std::mt19937 rng;
  if (opt.random_seed) {
    std::random_device rd;
    opt.seed = rd();
    rng.seed(opt.seed);
    std::cout << "Random seed: " << opt.seed << "\n";
  } else {
    rng.seed(opt.seed);
  }

  std::cout << "Grid " << opt.nx << " x " << opt.ny << " x " << opt.nz
            << "  h=" << opt.h_mm << " mm\n";
  std::cout << "Config: " << opt.config_path << "\n";

  try {
    ChemicalEnvironment env(opt.nx, opt.ny, opt.nz, opt.h_mm);
    env.load_from_config(opt.config_path, opt.diffusivity_scale, 1.0);
    env.allocate_channels_from_config();
    if (opt.algorithm_choice >= 0)
      env.set_diffusion_algorithm(
          static_cast<DiffusionAlgorithm>(opt.algorithm_choice));

    const double tick_dt = opt.tick_dt_minutes > 0.0
                               ? opt.tick_dt_minutes
                               : env.tick_interval_minutes();
    const double d_tnf = env.registry().diffusivity(TNF);
    const double spread_mm = std::sqrt(6.0 * d_tnf * tick_dt);

    std::cout << "Diffusion algorithm: " << env.diffusion_algorithm_label()
              << "\n";
    std::cout << "Diffusivity scale (Q_demo): " << opt.diffusivity_scale
              << "  D_eff(TNF) = " << d_tnf << " mm^2/min\n";
    std::cout << "Diffusion tick dt (min): " << tick_dt << "\n";
    std::cout << "Estimated RMS spread per tick ~ " << spread_mm << " mm (~"
              << (spread_mm / opt.h_mm) << " patches)\n";

    place_random_spikes(env, opt, rng);
    print_stats(env, opt.nx, opt.ny, opt.nz, "After spikes");

    if (opt.write_paraview) {
      std::error_code ec;
      std::filesystem::create_directories(opt.output_dir, ec);
      write_paraview_timestep(env, opt, 0);
    }

    for (int t = 0; t < opt.num_ticks; ++t) {
      env.run_diffusion_phase(tick_dt);
      env.merge_and_reset_secretion();
      env.update_chemotaxis_from_species(TGF);

      std::ostringstream label;
      label << "After tick " << (t + 1);
      print_stats(env, opt.nx, opt.ny, opt.nz, label.str().c_str());

      if (opt.write_paraview)
        write_paraview_timestep(env, opt, t + 1);
    }

    if (opt.write_paraview) {
      std::cout << "Open in ParaView: File → Open → " << opt.output_dir
                << "/chem_spikes_t*.vti"
                << " → group as time series.\n";
    }
    std::cout << "Done.\n";
  } catch (const std::exception &ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
