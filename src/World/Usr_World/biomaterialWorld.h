/*
 * biomaterialWorld.h
 *
 * File Contents: Contains declarations for the BMWorld class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#ifndef BMWORLD_H
#define BMWORLD_H

#include "../../Agent/Usr_Agents/Cell.h"
#include "../../ArrayChain/ArrayChain.h"
#include "../../Chemistry/chemical_environment.h"
#include "../../ECM/ECM.h"
#include "../../FieldVariable/Usr_FieldVariables/Chemical.h"
#include "../../common.h"
#include "../World.h"

#include <map>
#include <memory>
#include <new>
#include <stdlib.h>
#include <vector>

class Cell;
class Stem;
class Progen;
class NP;
class Collagen;
class Aggrecan;
class Hyaluronan;
class ECM;

using namespace std;

/*
 BMWORLD (BIOMATERIAL WORLD) CLASS DESCRIPTION: BMWorld is a derived class of
 the parent class World.
 *                                                The BMWorld class manages the
 model world.
 *                                                It is used to initialize
 cells, ECM, patches, and chemicals; to destroy
 *                                                agent ArrayChains; to execute
 each timestep of the model; to sprout agents;
 *                                                to count patches; and to
 output data.
 */
class BMWorld : public World {
public:
  /*
   * Description:	BMWorld constructor.
   *
   * Return: void
   *
   * Parameters: width    -- Width (x dimension) of the world in millimeters
   *             length   -- Length (y dimension) of the world in millimeters
   *             height   -- Height (z dimension) of the world in millimeters
   *             plength  -- Length of each patch (grid point) in millimeters
   */
  BMWorld(double width = 5,     // mm
          double length = 4,    // mm
          double height = 3,    // mm
          double plength = 0.01 // mm (10 um)
  );

  /*
   * Description:	BMWorld destructor.
   *
   * Return: void
   * Parameters: void
   * NOTE: Function is called implicitly: chonds.~ArrayChain()
   */
  ~BMWorld();

  /*
   * Description:	Destructor function for the Cell ArrayChain
   *
   * Return: void
   *
   * Parameters: &agent  -- Reference to cell that will be destroyed
   */
  void destroyCell(Cell *&agent);

  /*
   * Description:	Assign a patch type to each patch within bounds of type
   *
   * Return: void
   * Parameters: void
   */
  void assignPatches(int type, int xmin, int xmax, int ymin, int ymax, int zmin,
                     int zmax);

  /*
   * Description:	Assign a patch type to each patch in the world in row
   * major index manner
   *
   * Return: void
   * Parameters: void
   */
  void initializePatches();

  /** Copy baseline_total_mass from chemical_environment.json into baselineChem.
   */
  void sync_baseline_chem_from_config();

  /** Set baseline p* on patches via ChemicalEnvironment (after facade is
   * wired). */
  void initializeChemBaseline();

  /*
   * Description:	Initializes all cells to their correct patches
   *
   * Return: void
   * Parameters: void
   */
  void initializeCells();

  /*
   * Description:	Initializes collagen, aggrecan and hyaluronan to their
   * correct patches
   *
   * Return: void
   * Parameters: void
   */
  void initializeECM();

  /*
   * Description:	Initializes all damaged patches. Sprouts platelets and
   * fragments ECM on the damaged patches.
   *
   * Return: void
   * Parameters: void
   */
  void initializeDamage();

  /*
   * Description:	Each call to go() simulates 30 minutes or 'real-world'
   * time of the biological model.
   *
   * Return: 0 on success
   * Parameters: void
   */
  int go();

  /*
   * Description:	Entry function for sprouting cells
   *              Selects the appropriate sprouting function to apply
   *
   * Return: void
   *
   * Parameters: num          -- Number of cells to sprout
   *             patchType    -- Type of patches of sprout on
   *             agentType    -- Type of agent to sprout
   *
   *             Physical boundaries of the sprouting area/volume:
   *             xmin         -- Left
   *             xmax         -- Right
   *             ymin         -- Top
   *             ymax         -- Bottom
   *             zmin         -- Near
   *             zmax         -- Far
   *
   *             bloodOrTiss  -- Pass in true if should be sprouted in blood
   * 							   Pass in false if
   * should be sprouted in tissue (default)
   */
  void sproutAgent(int num, int patchType, int agentType, int xmin, int xmax,
                   int ymin, int ymax, int zmin, int zmax);

  /*
   * Description:	Function for sprouting cells in a given area/volume
   *
   * Return: void
   *
   * Parameters: num          -- Number of cells to sprout
   *             patchType    -- Type of patches of sprout on
   *             agentType    -- Type of agent to sprout
   *
   *             Physical boundaries of the sprouting area/volume:
   *             xmin         -- Left
   *             xmax         -- Right
   *             ymin         -- Top
   *             ymax         -- Bottom
   *             zmin         --
   *             zmax         --
   *
   *             bloodOrTiss  -- Pass in true if should be sprouted in blood
   *                             Pass in false if should be sprouted in tissue
   * (default)
   */
  void sproutAgentInArea(int num, int patchType, int agentType, int xmin,
                         int xmax, int ymin, int ymax, int zmin, int zmax);

  /*
   * Description:	Function for sprouting cells in the whole world
   *
   * Return: void
   *
   * Parameters: num          -- Number of cells to sprout
   *             patchType    -- Type of patches of sprout on
   *             agentType    -- Type of agent to sprout
   *             bloodOrTiss  -- Pass in true if should be sprouted in blood
   *                             Pass in false if should be sprouted in tissue
   * (default)
   */
  void sproutAgentInWorld(int num, int patchType, int agentType);

  /*
   * Description:	Counts the number of patches of a given type
   *
   * Return: int  -- Number of patches of given patch type
   *
   * Parameters: int  -- Enumic value for the patch type to count
   */
  int countPatchType(int);

  /****************************************************************
   * HELPER SUBROUTINES                                           *
   ****************************************************************/

  /*
   * Description:	Converts a length in millimeters to the number of
   * patches
   *
   * Return: The number of patches
   *
   * Parameters: mm  -- Length in millimeters
   */
  int mmToPatch(double mm);

  /*
   * Description:	Converts hours and days into ticks
   *
   * Return: The number of ticks
   *
   * Parameters: hour  -- Number of hours
   *             day   -- Number of days
   */
  static int reportTick(int hour = 0, int day = 0);

  /*
   * Description:	Determines the number of minutes elapsed
   *
   * Return: The number of minutes elapsed
   *
   * Parameters: void
   */
  static double reportMinute();

  /*
   * Description:	Determines the number of hours elapsed
   *
   * Return: The number of hours elapsed
   *
   * Parameters: void
   */
  static double reportHour();

  /** Chemistry facade (null before constructor finishes chem init). */
  ChemicalEnvironment *chemical_environment();
  const ChemicalEnvironment *chemical_environment() const;

  /** Integrated cytokine masses (from environment when available). */
  float world_total_tnf() const;
  float world_total_tgf() const;
  float world_total_il1beta() const;

  /** Integrated O2 mass */
  float world_total_o2() const;

  float chem_concentration(SpeciesId species, int patch_index) const;
  void chem_add_secretion(SpeciesId species, int patch_index,
                          float delta) const;
  float chem_concentration_channel(int concentration_channel,
                                   int patch_index) const;
  float chemotaxis_at(int patch_index) const;

  /*
   * Description:	Determines the number of days elapsed
   *
   * Return: The number of days elapsed
   *
   * Parameters: void
   */
  static double reportDay();

  /*
   * Description:	Determines the number of neighbors with given patch type
   *
   * Return: The number of neighbors with given patch type
   *
   * Parameters: ix         -- x coordinate of current patch
   *             iy         -- y coordinate of current patch
   *             iz         -- z coordinate of current patch
   *             patchType  -- patch type to count
   */
  int countNeighborPatchType(int ix, int iy, int iz, int patchType);

  /*
   * Description:	Reads user input from config file (default: config.txt)
   * to initialize chemicals, wound, cells.
   *
   * Return: Returns 0 if function proceeded to completion, for testing. Can be
   * removed later.
   *
   * Parameters: void
   */
  int userInput();

  // #ifdef MODEL_SCAFFOLD
  /*
   * Description: Calculate Hydrogel initial Elastic modulus, Crosslink density,
   * pore size
   *
   * Return:
   * Parameters:
   */
  void initializeCaAlg();

  /*
   * Description: Update Ca-Alg Swelling Ratio at current tick
   *
   * Return:
   * Parameters:
   */
  void updateSwellingRatio();

  /*
   * Description:Update Ca-Alg Mass Loss (% of initial mass) at current tick
   *
   * Return:
   * Parameters:
   */
  void updateMassLoss();

#ifdef PEPTIDE_BM
  /*
   * Description: Update effective stiffness E from peptide stress-relaxation
   * parameters (E_0, E_inf, t).
   */
  void updateE();
#endif

  /*
   * Description:Degrade Ca-Alg Patches and replace with immature tissue patch
   * (no new ECM produced)
   *
   * Return:
   *
   * Parameters:        numOfPatches  -- number of CaAlg patches to "degrade"
   * and be replaced with patch type tissue
   */
  void degradeCaAlg(int numOfPatches);

  /*
   * Description: Print out extra info for debugging purposes
   *
   * Return:
   * Parameters:
   */
  void debugInfo();
  // #endif

  /****************************************************************
   * OUTPUT SUBROUTINES & VISUALIZATION                           *
   ****************************************************************/

  ///*
  // * Description:	Outputs cell counts and cytokine levels from the current
  // tick to the file "Output/Output_Biomarkers.csv".
  // *              Used for testing.
  // *
  // * Return: void
  // * Parameters: void
  // */
  // void outputWorld_csv();

  /*
   * Description:	Outputs all patch assignments (patch type, agent type,
   * ECM type) to files in output directory.
   *
   * Return: void
   * Parameters: void
   */
  void patchassign_csv();

  /****************************************************************
   * STATIC VARIABLES                                             *
   ****************************************************************/

  static unsigned seed;      // Used to generate random numbers
  static bool highTNFdamage; // Whether there is high TNF damage (which results
                             // in ECM fragmentation)
  static float patchpermm; // The number of patches per millimeter in the world
  static float liveCells;
  static float deadCells;
  static float deletedCells;
  static float prevCells;
  static int initialCaAlg; // The number of initial tissue patches
  
  static float E; // Effective stiffness

  static float initialO2; // Initial concentration of oxygen (umol/L)

#ifdef PEPTIDE_BM
  static float E_0;   // Initial elastic modulus (peptide-conjugated BM)
  static float E_inf; // Equilibrium modulus (peptide-conjugated BM)
  static float t;     // Stress relaxation time (seconds)
#endif

  // Crosslinked Ca-Alg Hydrogel Parameters
  //  Ca Crosslinker:
  static float Ca_Mw; // Molecular Weight
  // Alginate:
  static float Alg_Mn;        // Number Average Molecular Weight (g mol-1)
  static float totalVolumeML; //

  /* CALIBRATION Variables */
  static float thresholdTNFdamage;    // The threshold for TNF damage
  static float sproutingFrequency[6]; // The number of hours between agent
                                      // sprouting sessions
  static float sproutingAmount[14]; // Constants related to the number of agents
                                    // to sprout
  static float cytokineDecay[6];    // The decay rates of the cytokines
  static float
      halfLifes_static[6]; // The half lifes of the cytokines in minutes

  /* Calibration Variables */
  static float ElasticMod[7]; // Elastic Modulus of Ca-Alg Hydrogel
  static float XLDensity[2];  // Crosslink density of Ca-Alg Hydrogel
  static float SwellRatio[5]; // Swell Ratio of Ca-Alg Hydrogel
  static float MassLoss[4];   // Mass loss of Ca-Alg Hydrogel
  static float PoreSize[2];   // PoreSize of Ca-Alg Hydrogel

  /****************************************************************
   * CONSTANT VARIABLES                                           *
   ****************************************************************/
  double patchlength; // The length of each patch

  /** World-level cytokine totals (grid owned by chemical_environment_). */
  Chemical WorldChem;

  float pXL;       // Crosslink Density (mmol/mL)
  float Q;         // Swelling Ratio (%)
  float w;         // Mass Loss (%)
  float poreWidth; // Pore Size (um)
  
#ifdef PEPTIDE_BM
  string peptide; // Type of peptide used for biomaterial conjugation
#endif

  int typesOfChem; // The number of different chemicals there are in the world
  vector<float> baselineChem; // Initial amount of each chemical in the world
  ECM *worldECM;              // Pointer to array of ECM
  vector<int>
      initialCells; // Initial amount of each cell type (agent) in the world

  ArrayChain<Cell *> cells; // ArrayChain to manage all cell data
  vector<Cell *>
      *localNewCells[MAX_NUM_THREADS]; // Vector of pointers to local lists of
                                       // cell pointers to add to global list
  vector<int> initHAcenters;       // Vector of patches which can be centers for
                                   // sprouting original hyaluronan
  unsigned seeds[MAX_NUM_THREADS]; // Seeds used to generate random numbers for
                                   // each thread

  vector<float> tgfLine; // Vector to store TGF values along an x-face line from
                         // boundary to center of ABM grid

  int lineY;
  int lineZ;

  float Alg_v, Alg_wv; // Volume (mL) and final concentration (% w/v) of Alg in
                       // Ca-Alg hydrogel
  float Ca_v, Ca_wv;   // Volume (mL) and final concentration (% w/v) of Ca 3400
  float highMW_alg, lowMW_alg; // ratio components of high and low MW kDa in the
                               // alginate hydrogel

  /** Owns patch grids, registry, and per-tick diffusion. */
  std::unique_ptr<ChemicalEnvironment> chemical_environment_;

  /** Fallback tick length (minutes) if environment not initialized. */
  static constexpr double kTickIntervalMinutes = 30.0;

  double tick_interval_minutes() const;

protected:
  // --- output file hooks ---
  char* get_output_filename() override;
  void write_csv_header(std::ofstream &file) override;
  void write_data_row(std::ofstream &file,
                      std::map<std::string, int> &agent_counts,
                      std::map<std::string, float> &env_counts) override;

  // --- auxiliary output hooks ---
  void write_auxiliary_header() override;
  void write_auxiliary_outputs() override;

  // --- agent counting hooks ---
  std::vector<std::string> get_agent_type_names() override;
  void count_agent_types(std::map<std::string, int> &agent_counts) override;

  // --- agent population hooks ---
  int get_total_agent_count() override;

  // --- environment element counting hooks (e.g. ecm) ---
  std::vector<std::string> get_env_type_names() override;
  void count_env(std::map<std::string, float> &env_counts) override;

private:
  /****************************************************************
   * MAJOR SECTION SUBROUTINES - begin                            *
   ****************************************************************/

  /*
   * Description: (Stage 0)	Sprout agents on various patches.
   *
   * Return: void
   * Parameters: hours  -- Current hour in model execution
   */
  // void seedCells(float hours);

  /*
   * Description: (Stage 1) Diffuse all registry species over one tick.
   *
   * Clears d*, reads p*, writes diffusion increment into d* (cells add
   * secretion later). PDE numerics live in ChemicalEnvironment /
   * diffusion3d_core ? not in this file.
   *
   * Return: void
   * Parameters: void
   */
  void diffuseCytokines();

  /*
   * Description:	Helper function for cell execution. Execute all alive
   * chondrocytes.
   *
   * Return: void
   *
   * Parameters: void
   */
  void inline runCells();

  /*
   * Description:	(Stage 2)	Execute cell functions for all living
   * cells.
   *
   * Return: void
   *
   * Parameters: void
   */
  void executeCells();

  /*
   * Description:	(Stage 3)	Execute ECM functions.
   *
   * Return: void
   *
   * Parameters: void
   */
  void executeECMs();

  /*
   * Description:	(Stage 3)	Fragment ECM proteins if necessary.
   *
   * Return: void
   *
   * Parameters: void
   */
  void requestECMfragments();

  /*
   * Description: (Stage 4a) Commit per-tick chemical changes to patch
   * concentrations.
   *
   * For each diffusing species: p* += d* (diffusion increment + cell
   * secretion), then d* = 0. Called after executeCells() so d* holds both
   * contributions when merged.
   *
   * Return: void
   * Parameters: void
   */
  void updateChem();

  /*
   * Description: Same merge as updateChem(); CPU implementation used each tick.
   *
   * Return: void
   * Parameters: void
   */
  void updateChemCPU();

  /*
   * Description:	Helper function for ECM updates. Execute updates for ALL
   * ECM managers.
   *
   * Return: void
   *
   * Parameters: void
   */
  void inline executeAllECMUpdates();

  /*
   * Description:	Helper function for ECM updates. Execute request resets
   * for ALL ECM managers.
   *
   * Return: void
   *
   * Parameters: void
   */
  void inline executeAllECMResetRequests();

  /*
   * Description:	(Stage 4c)	Update ECM managers to reflect next
   * tick's states
   *
   * Return: void
   *
   * Parameters: void
   */
  void updateECMManagers();

  /*
   * Description:	(Stage 4b)	Update cells to reflect next tick's
   * states
   *
   * Return: void
   *
   * Parameters: void
   */
  void updateCells();

  /*
   * Description:	Update cells to reflect next tick's states.
   *              This is called instead of updateCells() during setup.
   *
   * Return: void
   *
   * Parameters: void
   */
  void updateCellsInitial();

  // --- output function hooks ---
  // viability and differentiation are internal calculations
  // used only in write_data_row ? not exposed as hooks
  float calculate_viability();
  float calculate_pct_differentiated(std::map<std::string, int> &agent_counts);
};

#ifdef PEPTIDE_BM
/* Experimental values for peptide biomaterial */
struct peptideCondition {
    float E_init; // initial modulus, represents elastic portion of model
    float E_eq; // 'equilibrium' modulus, represents viscous portion of model
    float t_stress; // stress relaxation time (t = t_half/ln(2))
};
extern peptideCondition MAL, CHAD, hA5G26, IKVAV; // names of the peptide structs
//struct peptideCondition MAL, CHAD, hA5G26, IKVAV; // names of the peptide structs
#endif //PEPTIDE_BM

#endif /* BMWorld_H */
