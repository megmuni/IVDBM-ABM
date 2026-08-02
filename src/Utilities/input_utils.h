/*
 * input_utils.h
 *
 * File Contents: Contains functions to manage user input for model properties
 *
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#ifndef INPUT_UTILS_H_
#define INPUT_UTILS_H_

#pragma once

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <sys/stat.h>

using namespace std;
namespace util {

/*************************************************************************
 * MODEL OPTIONS                                                         *
 *************************************************************************/
int numTicks; // Desired number of ticks for the simulation (1 tick = 30 min)
float patchWidth;  // Desired width of each patch in millimeters
float worldXwidth; // Desired x-dimension width of the world in millimeters
float worldYwidth; // Desired y-dimension width of the world in millimeters
float worldZwidth; // Desired z-dimension width of the world in millimeters
/** Path to simulation_config.json (copy from simulation_config.template.json).
 */
char simulationConfigFile[200];
char outputDir[200];      // Base directory for run artifacts
char outputFileName[200]; // Primary biomarker CSV path
bool paraviewExportEnabled = false;

inline const char *getOutputDir() { return outputDir; }

inline bool paraviewEnabled() { return paraviewExportEnabled; }

inline void makeOutputPath(char *dest, size_t dest_size, const char *filename) {
  snprintf(dest, dest_size, "%s/%s", outputDir, filename);
}

inline void copyArgToBuffer(char *dest, size_t dest_size, const char *arg) {
  strncpy(dest, arg, dest_size - 1);
  dest[dest_size - 1] = '\0';
}

inline void mkdir_p(const char *path) {
  if (path == nullptr || path[0] == '\0')
    return;

  char buf[512];
  strncpy(buf, path, sizeof(buf) - 1);
  buf[sizeof(buf) - 1] = '\0';

  for (char *p = buf + 1; *p != '\0'; ++p) {
    if (*p == '/') {
      *p = '\0';
      mkdir(buf, 0755);
      *p = '/';
    }
  }
  mkdir(buf, 0755);
}

inline void ensureOutputDir() { mkdir_p(outputDir); }

inline void ensureOutputSubpath(const char *subpath) {
  char path[512];
  makeOutputPath(path, sizeof(path), subpath);
  mkdir_p(path);
}

inline void initOutputDirFromEnv() {
  const char *env = getenv("IVDBM_OUTPUT_DIR");
  if (env == nullptr || env[0] == '\0')
    env = getenv("OUTPUT_DIR");
  if (env != nullptr && env[0] != '\0')
    strncpy(outputDir, env, sizeof(outputDir) - 1);
  else
    strncpy(outputDir, "output", sizeof(outputDir) - 1);
  outputDir[sizeof(outputDir) - 1] = '\0';
}

inline void finalizeOutputPaths() {
  if (outputFileName[0] == '\0')
    snprintf(outputFileName, sizeof(outputFileName), "%s/Output_Biomarkers.csv",
             outputDir);
}

/*
 * Description: Path to the simulation JSON config (world_init + chemistry +
 * biology sections).
 */
inline const char *getSimulationConfigPath() { return simulationConfigFile; }

/*
 * Description:	Function for getting the number of ticks that the user inputted
 *
 * Return: The number of ticks that the user inputted
 * Parameters: void
 */
int getNumTicks() { return numTicks; }

/*
 * Description:	Function for getting the patch width that the user inputted
 *
 * Return: The patch width that the user inputted
 * Parameters: void
 */
float getPatchWidth() { return patchWidth; }

/*
 * Description:	Function for getting the world x,y,z -widths that the user
 * inputted
 *
 * Return: The width of the world in the x, y, z dimensions that the user
 * inputted Parameters: void
 */
float getWorldXWidth() { return worldXwidth; }
float getWorldYWidth() { return worldYwidth; }
float getWorldZWidth() { return worldZwidth; }

/*
 * Description:	Function for turning command line arguments into model options
 *
 * Return: void
 * Parameters: argc  -- Command line arguments count
 *             argv  -- Command line arguments
 */
void processOptions(int argc, char **argv) {

  outputFileName[0] = '\0';
  initOutputDirFromEnv();

/* ------------------------- Setting default options ------------------------ */
#ifdef MODEL_SCAFFOLD
  numTicks = 432;
#else
  numTicks = 500;
#endif

#ifdef MODEL_VOCALFOLD
  patchWidth = 0.010; // (mm)
  worldXwidth = 24.9; // (mm)
  worldYwidth = 1.6;  // (mm)
#ifdef MODEL_3D
  worldZwidth = 17.4; // (mm)
#else
  worldZwidth = patchWidth;
#endif
#elif defined(MODEL_SCAFFOLD)
  // Default 0.03 mL Scaffold with same cell seeding density as input
  patchWidth = 0.01; // (mm)
  worldXwidth = 3.1; // (mm)
  worldYwidth = 3.1; // (mm)
#ifdef MODEL_3D
  worldZwidth = 3.1; // (mm)
#else
  worldZwidth = patchWidth; // (mm)
#endif
#endif // MODEL_VOCALFOLD

  strcpy(simulationConfigFile, "configFiles/simulation_config.json");

  if (argc == 1) {
    finalizeOutputPaths();
    return;
  }

  // Get options from command line arguments
  for (int i = 1; i < argc; i++) {
    char *option_string = argv[i];
    if (!strcmp(option_string, "--numticks")) {
      numTicks = atoi(argv[++i]);
      //} else if (!strcmp(option_string, "--patchwidth")) {    // patch width
      // fixed at 0.01mm
      // maintain 1 cell max per patch occupancy
      // patchWidth = atof(argv[++i]);
    } else if (!strcmp(option_string, "--wxw")) {
      worldXwidth = atof(argv[++i]);
    } else if (!strcmp(option_string, "--wyw")) {
      worldYwidth = atof(argv[++i]);
    } else if (!strcmp(option_string, "--wzw")) {
#ifdef MODEL_3D
      worldZwidth = atof(argv[++i]);
#ifdef PDE_DIFFUSE
      if (worldZwidth < 0.06) { // minimum z dimension for 3D PDE diffuseChem
        cerr << " Error: 3D wzw must be greater than 0.06mm" << endl;
        exit(1);
      }
#endif
#else
      cerr << "Error: 3D functionalities undefined. Enter 2D World dimensions. "
           << endl;
      exit(1);
#endif

    } else if (!strcmp(option_string, "--outputfile")) {
      copyArgToBuffer(outputFileName, sizeof(outputFileName), argv[++i]);
    } else if (!strcmp(option_string, "--output-dir")) {
      copyArgToBuffer(outputDir, sizeof(outputDir), argv[++i]);
    } else if (!strcmp(option_string, "--config")) {
      copyArgToBuffer(simulationConfigFile, sizeof(simulationConfigFile),
                      argv[++i]);
    } else if (!strcmp(option_string, "--paraview")) {
      paraviewExportEnabled = true;
    } else if (!strcmp(option_string, "--help")) {
      cout << "Options: " << endl;
      cout << "   --numticks:      Number of ticks" << endl;
      // cout << "   --patchwidth:    Patch width    (mm)" << endl;
      cout << "   --wxw:           World width    (mm)" << endl;
      cout << "   --wyw:           World length   (mm)" << endl;
      cout << "   --wzw:           World height   (mm)" << endl;
      cout << "   --outputfile:    path to biomarker CSV (default: "
              "<output-dir>/Output_Biomarkers.csv)"
           << endl;
      cout << "   --output-dir:    directory for run output (default: output, "
              "or $IVDBM_OUTPUT_DIR)"
           << endl;
      cout << "   --config:        path to simulation_config.json (world_init + "
              "chemistry + biology)"
           << endl;
      cout << "   --paraview:      export patch/chemistry .vti time series "
              "under <output-dir>/paraview/"
           << endl;
      cout << " Usage: " << endl;
      cout << "      ./bin/testRun --wxw 3.1 --wyw 3.1 --wzw 3.1 --numticks 48"
           << endl;
      exit(-1);
    } else {
      cerr << "Error: Invalid option: " << option_string << endl;
      exit(-1);
    }
  }

  finalizeOutputPaths();
}

inline std::string jsonEscape(const char *s) {
  if (s == nullptr)
    return std::string();

  std::string out;
  out.reserve(strlen(s) + 8);
  for (const char *p = s; *p != '\0'; ++p) {
    switch (*p) {
    case '"':
      out += "\\\"";
      break;
    case '\\':
      out += "\\\\";
      break;
    case '\n':
      out += "\\n";
      break;
    case '\r':
      out += "\\r";
      break;
    case '\t':
      out += "\\t";
      break;
    default:
      if (static_cast<unsigned char>(*p) < 0x20) {
        char buf[8];
        snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned char>(*p));
        out += buf;
      } else {
        out += *p;
      }
    }
  }
  return out;
}

inline void writeJsonStringField(std::ofstream &out, const char *key,
                                 const char *value, bool last = false) {
  out << "    \"" << key << "\": ";
  if (value != nullptr && value[0] != '\0')
    out << "\"" << jsonEscape(value) << "\"";
  else
    out << "null";
  if (!last)
    out << ",";
  out << "\n";
}

/*
 * Description: Write resolved run parameters to JSON when IVDBM_RUN_PARAMS_JSON
 *              is set (Slurm jobs via scripts/testrun.sbatch).
 */
inline void writeRunParamsJson(int argc, char **argv) {
  const char *path = getenv("IVDBM_RUN_PARAMS_JSON");
  if (path == nullptr || path[0] == '\0')
    return;

  std::ofstream out(path);
  if (!out) {
    fprintf(stderr, "Warning: cannot write run params to %s\n", path);
    return;
  }

  out << "{\n";
  out << "  \"slurm\": {\n";
  writeJsonStringField(out, "job_id", getenv("SLURM_JOB_ID"));
  writeJsonStringField(out, "array_task_id", getenv("SLURM_ARRAY_TASK_ID"));
  writeJsonStringField(out, "job_name", getenv("SLURM_JOB_NAME"));
  writeJsonStringField(out, "cpus_per_task", getenv("SLURM_CPUS_PER_TASK"));
  writeJsonStringField(out, "mem_per_node", getenv("SLURM_MEM_PER_NODE"));
  writeJsonStringField(out, "submit_dir", getenv("SLURM_SUBMIT_DIR"));
  writeJsonStringField(out, "partition", getenv("SLURM_JOB_PARTITION"));
  writeJsonStringField(out, "qos", getenv("SLURM_JOB_QOS"));
  writeJsonStringField(out, "nodelist", getenv("SLURM_JOB_NODELIST"));
  writeJsonStringField(out, "submit_host", getenv("SLURM_SUBMIT_HOST"));
  writeJsonStringField(out, "gpus_on_node", getenv("SLURM_GPUS_ON_NODE"));
  writeJsonStringField(out, "account", getenv("SLURM_JOB_ACCOUNT"));
  writeJsonStringField(out, "profile", getenv("IVDBM_PROFILE"), true);
  out << "  },\n";

  out << "  \"simulation\": {\n";
  out << "    \"num_ticks\": " << numTicks << ",\n";
  out << "    \"patch_width_mm\": " << patchWidth << ",\n";
  out << "    \"world_x_width_mm\": " << worldXwidth << ",\n";
  out << "    \"world_y_width_mm\": " << worldYwidth << ",\n";
  out << "    \"world_z_width_mm\": " << worldZwidth << ",\n";
  writeJsonStringField(out, "output_dir", outputDir);
  writeJsonStringField(out, "output_file", outputFileName);
  writeJsonStringField(out, "simulation_config", simulationConfigFile);
  out << "    \"paraview_enabled\": "
      << (paraviewExportEnabled ? "true" : "false") << "\n";
  out << "  },\n";

  out << "  \"command\": {\n";
  out << "    \"argv\": [";
  for (int i = 0; i < argc; ++i) {
    if (i > 0)
      out << ", ";
    out << "\"" << jsonEscape(argv[i]) << "\"";
  }
  out << "]\n";
  out << "  }\n";
  out << "}\n";
}

/*
 * Description: Write per-tick simulation timing when IVDBM_RUN_TIMING_JSON is
 *              set (merged into run_params.json after the job by sbatch).
 */
inline void writeRunTimingJson(long total_tick_ms, int num_ticks,
                               double setup_wall_seconds) {
  const char *path = getenv("IVDBM_RUN_TIMING_JSON");
  if (path == nullptr || path[0] == '\0')
    return;

  std::ofstream out(path);
  if (!out) {
    fprintf(stderr, "Warning: cannot write run timing to %s\n", path);
    return;
  }

  const double avg_ms =
      num_ticks > 0 ? static_cast<double>(total_tick_ms) / num_ticks : 0.0;

  out << "{\n";
  out << "  \"num_ticks_completed\": " << num_ticks << ",\n";
  out << "  \"setup_wall_seconds\": " << setup_wall_seconds << ",\n";
  out << "  \"tick_execution_ms_total\": " << total_tick_ms << ",\n";
  out << "  \"tick_execution_ms_avg\": " << avg_ms << "\n";
  out << "}\n";
}

/*
 * Description: Function for displaying model options
 *
 * Return: void
 * Parameters: void
 */
void printOptions() {
  cout << "ABM Parameters:" << endl;
  cout << "	numTicks:	" << numTicks << endl;
  cout << "	patchWidth:	" << patchWidth << " mm" << endl;
  cout << "	worldXwidth:	" << worldXwidth << " mm" << endl;
  cout << "	worldYwidth:	" << worldYwidth << " mm" << endl;
  cout << "	worldZwidth:	" << worldZwidth << " mm" << endl;
  cout << "	outputDir:	" << outputDir << endl;
  cout << "	outputFileName:	" << outputFileName << endl;
  cout << "	simulationConfig:	" << simulationConfigFile << endl;
  cout << "	paraviewEnabled:	"
       << (paraviewExportEnabled ? "true" : "false") << endl;
}

} // namespace util

#endif /* INPUT_UTILS_H_ */