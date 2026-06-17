/*
 * biomaterialWorld.cpp
 *
 * File contents: Contains the BMWorld class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <tgmath.h>
#include <time.h>
#include <vector>
#define PI 3.14159
#include "../../Utilities/error_utils.h"
#include "../../Utilities/input_utils.h"
#include "../../Utilities/timer.h"
#include "../../enums.h"
#include "biomaterialWorld.h"
#include <omp.h>
#include <sstream>
#include <string>

using namespace std;

/* ------------------------------------------------------------------------------------
 */
/*                            STATIC VARIABLES INITIALIZATIONS */
/* ------------------------------------------------------------------------------------
 */
unsigned BMWorld::seed = 27000; // initial number of cells
bool BMWorld::highTNFdamage = false;
float BMWorld::patchpermm = 0;
float BMWorld::liveCells = 0;
float BMWorld::deadCells = 0;
float BMWorld::deletedCells = 0;
float BMWorld::prevCells = 0;
#ifdef MODEL_SCAFFOLD
int BMWorld::initialCaAlg = 0;
float BMWorld::E = 0;     // Effective stiffness
float G;                  // Elastic Modulus (kPa)
float pXL;                // Crosslink Density (mmol/mL = M)
float Alg_Mn;             // molecular weight of alginate (kDa)
int highMW_alg;           // ratio component of high-MW alginate
int lowMW_alg;            // ratio component of low-MW alginate
float Q;                  // Swelling Ratio
float w = 0;              // Mass Loss (%)
float poreWidth = 200.00; // (um)
#endif

#ifdef PEPTIDE_BM
	string peptide; // Type of peptide used for biomaterial conjugation
	float BMWorld::E_0 = 0;
	float BMWorld::E_inf = 0;
	float BMWorld::t = 0;
#endif

#ifdef MODEL_SCAFFOLD
float BMWorld::Ca_Mw = 3400; // Ca Molecular Weight (Mw ≈ 3,400 = g/mol)
float BMWorld::Alg_Mn =
    1500; // Average molecular weight (Mw = 1 kDa = 1000 g/mol)
float BMWorld::totalVolumeML;
// float BMWorld::Alg_Mn = 90;
// float BMWorld::Alg_Mn = 200; //143;
#endif

float BMWorld::thresholdTNFdamage = 10.0; // ng //unused in IVDBM-ABM
float BMWorld::cytokineDecay[6] = {0.2, 0.2, 0.2, 0.2, 0.2, 0.5};  // 0.2, 0.2,
float BMWorld::halfLifes_static[6] = {33.6, 2.7, 46, 103, 24, 60}; // 13, 13,

#ifdef MODEL_SCAFFOLD
float BMWorld::ElasticMod[7] = {125, 58, 971, 1.037, 756, 0.516, 0.165};
float BMWorld::XLDensity[2] = {
    2.3, 10.1}; // IN IVDBM-ABM (stem cell version) THESE ARE NOT USED ANYWHERE
float BMWorld::SwellRatio[5] = {
    72.478, 0.131, 22.034, 3.284,
    35.752}; // float BMWorld::SwellRatio[5] = {0.4, 0.4, 3, 7.9, 1400};
float BMWorld::MassLoss[4] = {
    0.234, 7.785, 0.15,
    1.36}; // float BMWorld::MassLoss[4] = {17.6, 0.9, 60, 5.3};
float BMWorld::PoreSize[2] = {
    1769.8, 258.5}; // float BMWorld::PoreSize[3] = {345.2, 309.9, 138.1};
#endif

#ifdef PEPTIDE_BM
	/* Experimental values for peptide biomaterial */
	peptideCondition MAL = { 2.0031, 3.8757, 44.00 };
	peptideCondition CHAD = { 3.0180, 5.0709, 37.87 };
	peptideCondition hA5G26 = { 2.9628, 4.6913, 35.35 };
	peptideCondition IKVAV = { 2.9067, 5.4795, 42.56 };
#endif

BMWorld::BMWorld(double length, double width, double height, double plength) {
  // Generate random seeds:
  for (int i = 0; i < NUM_THREAD; i++)
    seeds[i] = 25234 + 17 * i;

  // Allocate memory for local lists of cell pointers to add:
  for (int i = 0; i < MAX_NUM_THREADS; i++)
    localNewCells[i] = new vector<Cell *>;

  /* --------------------------------------------------------------------------
   */
  /*                                 GRID SETUP */
  /* --------------------------------------------------------------------------
   */
  this->patchlength = plength;

  // Number of patches in x,y,z dimensions:
  int nx = width / patchlength;
  int ny = length / patchlength;
  int nz = height / patchlength;
  World::setupGrid(nx,     // number of grid points (patches) in x dimension
                   ny,     // number of grid points (patches) in y dimension
                   nz,     // number of grid points (patches) in z dimension
                   0.0,    // min coordinates in x
                   width,  // max corodinates in x
                   0.0,    // min coordinates in y
                   length, // max coordinates in y
                   0.0,    // min coordinates in z
                   height  // max coordinates in z
  );

  // Read input parameters (chem baseline, wound dimensions, initial cells) from
  // config file
  int temp = BMWorld::userInput();
  cout << "length, width, height: " << length << " mm, " << width << " mm, "
       << height << " mm" << endl;
  cout << "Number of patches: nx, ny, nz: " << nx << ", " << ny << ", " << nz
       << " " << endl;

  // Allocate and initialize Patches/ECM
  if (util::ABMerror(!(this->worldPatch = new Patch[(nx) * (ny) * (nz)]),
                     "Patch mem alloc error!", __FILE__, __LINE__))
    exit(1);
  if (util::ABMerror(!(this->worldECM = new ECM[(nx) * (ny) * (nz)]),
                     "ECM mem alloc error!", __FILE__, __LINE__))
    exit(1);
  cout << "worldPatch size: " << (nx) * (ny) * (nz)
       << " (Number of world patches) " << endl;
  cout << "worldECM size: " << (nx) * (ny) * (nz) << " (Number of ECM patches) "
       << endl;

  /* Try initializing Patches and ECMs with the threads that will access them
   * later since the default allocation policy on Linux platforms is
   * first-touch. This is a best-effort implementation, since we cannot
   * guarantee size of data accessed per thread to be an integer multiple of
   * page size. */
  std::cout << std::fixed;
  std::cout << std::setprecision(3);
  // cout << "	allocating ECM Managers (also Patches) with best-effort first
  // touch" << endl;

  for (int iz = 0; iz < nz; iz++) {
#pragma omp parallel for
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        this->worldPatch[in] = Patch(ix, iy, iz, in);
        this->worldECM[in] = ECM(ix, iy, iz, in);
      }
    }
  }

  // Define Class static variables and pointers
  BMWorld::patchpermm = nx / width;
  Agent::nx = this->nx;
  Agent::ny = this->ny;
  Agent::nz = this->nz;
  Agent::agentWorldPtr = this;
  Agent::agentPatchPtr = this->worldPatch;
  Agent::agentECMPtr = this->worldECM;
  ECM::ECMWorldPtr = this;
  ECM::ECMPatchPtr = this->worldPatch;

  /* ----------------------- INITIALIZATION SUBROUTINES -----------------------
   */

  /* Define initial attributes of patches, damage, ECM, chem and cells based on
   * user defined values (in config file) and traits of native tissue */
  this->initializePatches();
  this->initializeECM();
  this->initializeCells();
#ifdef MODEL_SCAFFOLD
  this->initializeCaAlg();
#endif
  /* Chemistry controller. */
  {
    const double swelling_Q =
        (this->Q > 0.0f) ? static_cast<double>(this->Q) : 1.0;
    chemical_environment_.reset(
        new ChemicalEnvironment(nx, ny, nz, this->patchlength));
    chemical_environment_->load_from_config(
        util::getChemicalEnvironmentConfigPath(), swelling_Q);
    chemical_environment_->allocate_channels_from_config();
    this->sync_baseline_chem_from_config();
    this->initializeChemBaseline();
  }
  this->initializeDamage();

  /* Calling update functions to synchronize read and write portion of the
   * attributes */
  // this->updateChemCPU();
  this->updateCellsInitial(); // Add cells to list before removal and updates
  this->updateECMManagers();
  this->updatePatches();
  // cout << "setupGrid complete" << endl;
  cout << "-------------------------------------------" << endl;
}

BMWorld::~BMWorld() {
  for (int i = 0; i < MAX_NUM_THREADS; i++)
    delete localNewCells[i];
  cerr << " removing dead cells" << endl;

  int cellsSize = cells.size();

#pragma omp parallel for
  for (int i = 0; i < cellsSize; i++) {
#ifdef _OMP
    int tid = omp_get_thread_num();
#else
    int tid = DEFAULT_TID;
#endif

    Cell *cell = cells.getDataAt(i);
    if (!cell)
      continue;
    cells.deleteData(i, tid);
    delete cell;
  }

  if (worldPatch != NULL)
    delete[] worldPatch;
  if (worldECM != NULL)
    delete[] worldECM;
  cout << "BMWorld has been successfully destructed." << endl;
}

void destroyCell(Cell *&agent) {
  if (agent) {
    free(agent);
    agent = NULL;
  }
}

void BMWorld::assignPatches(int type, int xmin, int xmax, int ymin, int ymax,
                            int zmin, int zmax) {
  // Assign patches within bounds of type 'type'
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        switch (type) {
        case damage:
          this->worldPatch[in].type[read_t] = damage;
          this->worldPatch[in].type[write_t] = damage;
          this->worldPatch[in].color[read_t] = cdamage;
          this->worldPatch[in].color[write_t] = cdamage;
          this->worldPatch[in].dirty = true;
          break;
        case CaAlg:
          this->worldPatch[in].type[read_t] = CaAlg;
          this->worldPatch[in].type[write_t] = CaAlg;
          this->worldPatch[in].color[read_t] = cCaAlg;
          this->worldPatch[in].color[write_t] = cCaAlg;
          this->worldPatch[in].dirty = true;
          break;
        }
      }
    }
  }
}

void BMWorld::initializePatches() {
#ifdef MODEL_3D
  assignPatches(CaAlg, 0, nx, 0, ny, 0, nz);
#else
  assignPatches(CaAlg, 0, nx, 0, ny, 0, 0);
#endif

  tgfLine.resize(nx / 2 + 1, 0.0f);
  lineY = ny / 2;
  lineZ = nz / 2;

  // Assign values to initial:
  BMWorld::initialCaAlg = this->countPatchType(CaAlg);
  // cout << "Finished building Ca-Alg Hydrogel" << endl;
}

#ifdef MODEL_SCAFFOLD
void BMWorld::initializeCaAlg() {
  cout << "Begin Calculating Ca-Alg Properties..." << endl;

  /* ---------------------- Parameters of Ca-Alg Scaffold ---------------------
   */
  float Alg_ww = this->Alg_wv / (this->Alg_wv);

  /* ----- Calculate initial bulk mechanical properties of Ca-Alg Scaffold ----
   */
  /* p_XL: Crosslink Density (mmol/mL = M)
   * 		 Linear dependence of Shear modulus on cross-link concentration
   * for constant polymer concentration
   */
  this->pXL = this->pXL / 1000; // convert mM to M
  // this->pXL = 0.014;

  if (this->highMW_alg == 1 && this->lowMW_alg == 0) { // 'high' condition
    this->Alg_Mn = 1500;
  } else if (this->highMW_alg == 0 && this->lowMW_alg == 1) { // 'low' condition
    this->Alg_Mn = 95;
  } else { // 'mix' condition: calculates a weighted avg molecular weight
    float highMW_kDa = 1500;
    float lowMW_kDa = 50;
    this->Alg_Mn =
        (pow(this->highMW_alg * highMW_kDa, 2) +
         pow(this->lowMW_alg * lowMW_kDa, 2)) /
        ((this->highMW_alg * highMW_kDa) + (this->lowMW_alg * lowMW_kDa));
  }

  cout << "		Final Alginate concentration (%w/v): " << this->Alg_wv
       << endl;
  cout << "       Alginate Molecular Weight (kDa) = " << this->Alg_Mn << endl;
  cout << "       Calcium Crosslinking Density (mmol/mL = M) = " << this->pXL
       << endl;

/* Calculate Initial Elastic Modulus E (kPa)
 *  E = a (( b*TotalProtein(w/v) + c)* Alg(w/w) + d*TP(w/v)) + e*(f*Alg(w/v) +
 * g)*XL(w/w)
 *
 *       Follows rule of mixtures where stiffness of mixture is weight average
 * of components. Linear dependence of modulus on cross-link concentration for
 * constant polymer concentration
 */
#ifdef CALIBRATION
  this->E = -BMWorld::ElasticMod[0] +
            BMWorld::ElasticMod[1] * (Alg_wv)-BMWorld::ElasticMod[2] * (pXL) +
            BMWorld::ElasticMod[3] * (Alg_Mn) +
            BMWorld::ElasticMod[4] * (Alg_wv) * (pXL)-BMWorld::ElasticMod[5] *
                (Alg_wv) * (Alg_Mn)-BMWorld::ElasticMod[6] * (pXL) * (Alg_Mn);
#else
  this->E = -125 + 58 * (Alg_wv)-971 * (pXL) + 1.037 * (Alg_Mn) +
            756 * (Alg_wv * pXL) - 0.516 * (Alg_wv * Alg_Mn) -
            0.165 * (pXL * Alg_Mn);
#endif

#ifdef PEPTIDE_BM
	if (this->peptide.compare("MAL") == 0) {
		BMWorld::E_0 = MAL.E_init;
		BMWorld::E_inf = MAL.E_eq;
		BMWorld::t = MAL.t_stress;
	}
	else if (this->peptide.compare("CHAD") == 0) {
	  BMWorld::E_0 = CHAD.E_init;
		BMWorld::E_inf = CHAD.E_eq;
		BMWorld::t = CHAD.t_stress;
	}
	else if (this->peptide.compare("hA5G26") == 0) {
		BMWorld::E_0 = hA5G26.E_init;
		BMWorld::E_inf = hA5G26.E_eq;
		BMWorld::t = hA5G26.t_stress;
	}
	else if (this->peptide.compare("IKVAV") == 0) {
		BMWorld::E_0 = IKVAV.E_init;
		BMWorld::E_inf = IKVAV.E_eq;
		BMWorld::t = IKVAV.t_stress;
	}
	cout << "Peptide: " << this->peptide << ". Parameters being used are: " << BMWorld::E_0 << BMWorld::E_inf << BMWorld::t << endl;
	BMWorld::E = BMWorld::E_inf + (BMWorld::E_0 - BMWorld::E_inf) * exp(-(BMWorld::clock * 30 * 60) / BMWorld::t); // converts tick to seconds
#endif

  cout << "       Elastic Modulus (kPa) = " << this->E << endl;

/* Pore Size (um): poreWidth = -a * Alg_ww^2 + b * Alg_ww + c */
#ifdef CALIBRATION
  this->poreWidth = -BMWorld::PoreSize[0] * (pXL) + BMWorld::PoreSize[1];
#else
  this->poreWidth =
      -1769.84 * (pXL) +
      258.5; //-0.3113*pow(Alg_ww,2) + 1.5*Alg_ww + 50;   //this->poreWidth =
             //(-0.01)*345.2*pow(Alg_ww,2) + 309.9*Alg_ww + 138.1;
#endif
  cout << "     this->poreWidth = -" << BMWorld::PoreSize[0] << "*" << (pXL)
       << " + " << BMWorld::PoreSize[1] << endl;
  cout << "       Pore Width (um): " << this->poreWidth << endl;

/* Swelling Ratio:
 *       Swelling Ratio increase with Alg content and with time
 *       Important in retaining water, facilitating diffusion.
 */
#ifdef CALIBRATION
  this->Q = BMWorld::SwellRatio[0] -
            BMWorld::SwellRatio[1] * (this->reportDay()) -
            BMWorld::SwellRatio[2] * (this->Alg_wv) -
            BMWorld::SwellRatio[3] * (this->reportDay()) * (this->pXL) +
            BMWorld::SwellRatio[4] * (this->Alg_wv) * (this->pXL);
#else
  this->Q = 72.478 - 0.131 * (this->reportDay()) -
            22.034 * (Alg_wv)-3.284 * (this->reportDay()) * (pXL) +
            35.752 * (Alg_wv) * (pXL);
#endif
  // cout << "    this->Q ="<< BMWorld::SwellRatio[0]<< "-"<<
  // BMWorld::SwellRatio[1]<<"*"<<(this->reportDay())<<" - "<<
  // BMWorld::SwellRatio[2]<<"*"<<(Alg_wv)<<" -"<<
  // BMWorld::SwellRatio[3]<<"*"<<(this->reportDay())<<"*"<<(pXL) <<" + "<<
  // BMWorld::SwellRatio[4]<<"*"<<(Alg_wv)<<"*"<<(pXL) << endl;
  cout << "       Swelling Ratio: " << this->Q << endl;

/* Mass Loss (%)
 *       Instable hydrogel degrades, replaced with cell-synthesized ECM proteins
 *       Mass loss fraction (%) increases with time and Alg content
 */
#ifdef CALIBRATION
  this->w = 0; // this->w = BMWorld::MassLoss[0] + BMWorld::MassLoss[1]*(pXL) +
               // BMWorld::MassLoss[2]*(reportDay()) -
               // BMWorld::MassLoss[3]*(pXL)*(reportDay());
#else
  this->w = 0.234 + 7.785 * (pXL) + 0.15 * (reportDay()) -
            1.36 * (pXL) * (reportDay());
#endif
  if (w < 0)
    w = 0; // no negative mass loss
  // cout << " 	this->w ="<< BMWorld::MassLoss[0]<< " +"<<
  // BMWorld::MassLoss[1]<<"*"<<(pXL)<<" +"<<
  // BMWorld::MassLoss[2]<<"*"<<(reportDay())<<" - "<<
  // BMWorld::MassLoss[3]<<"*"<<(pXL)<<"*"<<(reportDay()) << endl;
  cout << " Mass Loss (%): " << this->w << endl;

  cout << "Finished calculating initial Ca-Alg properties" << endl;
}
#endif // MODEL_SCAFFOLD

void BMWorld::sync_baseline_chem_from_config() {
  if (!chemical_environment_) {
    if (util::ABMerror(1, "ChemicalEnvironment not initialized", __FILE__,
                       __LINE__))
      exit(1);
  }

  const ChemicalEnvironmentConfig &cfg = chemical_environment_->configuration();
  this->typesOfChem = cfg.channel_count;

  this->baselineChem.assign(4, 0.f);
  this->baselineChem[TNF] =
      chemical_environment_->baseline_total_mass_for("TNF");
  this->baselineChem[TGF] =
      chemical_environment_->baseline_total_mass_for("TGF");
  this->baselineChem[IL1beta] =
      chemical_environment_->baseline_total_mass_for("IL1beta");

  cout << "Chemical environment config: " << cfg.model << " (schema "
       << cfg.schema_version << ")" << endl;
  cout << "  tick_interval_minutes = " << cfg.tick_interval_minutes << endl;
  cout << "  channel_count = " << cfg.channel_count << endl;
  cout << "  diffusion algorithm: "
       << chemical_environment_->diffusion_algorithm_label() << endl;
  const DiffusionAlgorithm requested =
      chemical_environment_->diffusion_algorithm();
  if (std::string(chemical_environment_->diffusion_algorithm_label()) !=
      diffusion_algorithm_label(requested)) {
    cout << "  diffusion requested: "
         << diffusion_algorithm_label(requested) << endl;
  }
}

void BMWorld::initializeChemBaseline() {
  if (!chemical_environment_) {
    if (util::ABMerror(1, "ChemicalEnvironment not initialized", __FILE__,
                       __LINE__))
      exit(1);
  }

  chemical_environment_->clear_delta_channels();
  this->WorldChem.totalTNF = 0;
  this->WorldChem.totalTGF = 0;
  this->WorldChem.totalIL1beta = 0;

  if (this->baselineChem.size() != 4) {
    if (util::ABMerror(1, "Error initializing chemicals!!", __FILE__, __LINE__))
      exit(1);
    return;
  }

  const int countCaAlg = this->countPatchType(CaAlg);
  const float tnf0 = this->baselineChem[TNF] / countCaAlg;
  const float tgf0 = this->baselineChem[TGF] / countCaAlg;
  const float il10 = this->baselineChem[IL1beta] / countCaAlg;

  for (int iz = 0; iz < this->nz; iz++) {
#pragma omp parallel for
    for (int iy = 0; iy < this->ny; iy++) {
      for (int ix = 0; ix < this->nx; ix++) {
        const int in = ix + iy * nx + iz * nx * ny;
        if (this->worldPatch[in].type[read_t] == CaAlg) {
          chemical_environment_->set_concentration(in, TNF, tnf0);
          chemical_environment_->set_concentration(in, TGF, tgf0);
          chemical_environment_->set_concentration(in, IL1beta, il10);
        } else {
          chemical_environment_->set_concentration(in, TNF, 0.f);
          chemical_environment_->set_concentration(in, TGF, 0.f);
          chemical_environment_->set_concentration(in, IL1beta, 0.f);
        }
      }
    }
  }

  chemical_environment_->update_chemotaxis_from_species(TGF);
  chemical_environment_->recompute_world_totals();
  chemical_environment_->copy_totals_to(this->WorldChem);

  cout << "		Initial cytokine concentrations: totalTNF = "
       << this->WorldChem.totalTNF
       << ", totalTGF = " << this->WorldChem.totalTGF
       << ", totalIL1beta = " << this->WorldChem.totalIL1beta << endl;
}

void BMWorld::initializeCells() {
  // cout << "Begin Initializing Cells..." << endl;

  // Instantiate Cell list:
  cells = ArrayChain<Cell *>(DEFAULT_DATA_SMALL, 4, NULL,
                             NULL); // BMWorld::destroyChond);
  // cout << "Initialize Cells..." << endl;

  // If initial cell count not input, seed scaffold with stem cells at density
  // 10^6 cell/mL
  double cellDensity; // cells/mm^3
  double scaffoldVolume = (nx * ny * nz) * pow(this->patchlength, 3); // mm^3
  double hydrogelVolume = scaffoldVolume * pow(10, -3);               // mL

  if (this->initialCells[0] == 0) {
    double cellpermL = 1 * pow(10, 6);     // (Xuan Li) 1 mill/ml
    cellDensity = cellpermL * pow(10, -3); // cells/mm^3
  } else {
    double cellpermL = this->initialCells[0] / scaffoldVolume;
    cellDensity = cellpermL * pow(10, -3);
  }

  this->initialCells[0] = cellDensity * scaffoldVolume;
  int initialScaffoldCells = this->initialCells[0];

  cout << "		Scaffold Volume: " << scaffoldVolume << " mm^3 ("
       << hydrogelVolume << " mL)" << endl;
  cout << "		Cell Density: " << cellDensity << " cells/mm^3  ("
       << 1.0 * pow(10, 6) << " cells/mL)" << endl;
  cout << " 		Seeding " << initialScaffoldCells
       << " cells in scaffold " << endl;

  BMWorld::totalVolumeML = hydrogelVolume;

  // Sprout cell seeded hydrogel with mesenchymal stem cells at density 10^6
  // cells/mL
  sproutAgent(initialScaffoldCells, // Number of cells to sprout
              CaAlg,                // Type of patch to sprout on
              stem,                 // Type of agent to sprout
              // Physical Boundary:
              0,  //  -- left
              nx, //  -- right
              0,  //  -- top
              ny, //  -- bottom
              0,  //  -- front
              nz  //  -- rear
  );
  prevCells = this->initialCells[0];
  // cout << "Finished initializing cells" << endl;
}

void BMWorld::initializeECM() {

  /* --------------------------------------------------------------------------
   */
  /*                                  COLLAGEN */
  /* --------------------------------------------------------------------------
   */
  // cout << "Begin Initializing Collagen..." << endl;
  sproutAgentInArea(nx * ny * nz, // Number of agents to sprout
                    CaAlg,        // patch type
                    new_coll,     // new collagen agent type
                    0,  // lowest x-coordinate agent could be sprouted at
                    nx, // highest x-coordinate agent
                    0,  // lowest y-coordinate
                    ny, // highest y-coordinate
                    0,  // lowest z-coordinate
                    nz  // highest z-coordinate
  );                    // collagen in CaAlg Scaffold

  /* --------------------------------------------------------------------------
   */
  /*                                  AGGRECAN */
  /* --------------------------------------------------------------------------
   */
  // cout << "Begin Initializing Aggrecan..." << endl;
  sproutAgentInArea(nx * ny * nz, // Number of agents to sprout
                    CaAlg,        // patch type
                    new_agg,      // new aggrecan agent type
                    0,  // lowest x-coordinate agent could be sprouted at
                    nx, // highest x-coordinate agent
                    0,  // lowest y-coordinate
                    ny, // highest y-coordinate
                    0,  // lowest z-coordinate
                    nz  // highest z-coordinate
  );                    // aggrecan in CaAlg Scaffold
}

void BMWorld::initializeDamage() {
  // Sprout Damage throughout Scaffold and fragment ECM proteins
  for (int iz = 0; iz < nz; iz++) {
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        worldPatch[in].inDamzone = true;
        worldPatch[in].health[write_t] = 0;
        worldPatch[in].damage[write_t] = worldPatch[in].damage[read_t]++;
        worldPatch[in].dirty = true;

        if (worldECM[in].empty[write_t] == true)
          continue;
        this->worldECM[in].fragmentNCollagen();
        this->worldECM[in].fragmentNAggrecan();
        worldPatch[in].updatePatch();
      }
    }
  }
  // cout <<  this->countPatchType(damage) << " damage created" << endl;
}

/*
 * Each call to BMWorld::go() performs the following major steps:
 * 	(0) Cell seedings
 * 	(1) Chemical diffusion - diffuseCytokines(): PDE step; writes increment
 * to d* (2) Cell function       - cells may add secretion into d* (3) ECM
 * function (4) Attributes synchronization a) Update chemicals - updateChem():
 * p* += d*, clear d* b) Update cells c) Update ECM managers d) Update patches
 */
int BMWorld::go() {
  cout << "-------------------------------------------" << endl;

  Patch *tempPatchPtr;
  Agent *tempAgentPtr;
  double hours = this->reportHour();
  double days = this->reportDay();

// Profiling options defined in common.h
#ifdef PROFILE_MAJOR_STEPS
  struct timeval start, end;
  long elapsed_time; // in milliseconds
#endif

  // Increment Clock in ticks (1 tick = 30 min)
  BMWorld::clock++;
  cout << "tick: " << clock << " , hour: " << hours << " , day: " << days
       << endl;

#ifdef PROFILE_MAJOR_STEPS
  /* TIME_STAGE is a macro for timing a command/function and printing the timing
   * info (See Utilities/time.h) */

  /* --------------------------- CHEMICAL DIFFUSION ---------------------------
   */
  TIME_STAGE(this->diffuseCytokines(), "Chemical diffusion", "1");

  /* ------------------------------ CELL SEEDING ------------------------------
   */
  // TIME_STAGE(this->seedCells(hours), "Cell seeding", "0");

  /* ------------------------------ CELL FUNCTION -----------------------------
   */
  TIME_STAGE(this->executeCells(), "Cells function", "2");

  /* ------------------------------ ECM FUNCTION ------------------------------
   */
  TIME_STAGE(this->executeECMs(), "ECM function", "3(a)");
  TIME_STAGE(this->requestECMfragments(), "ECM fragmentation",
             "3(b)"); // Request fragment<ECM>

  /* ----------------------- ATTRIBUTES SYNCHRONIZATION -----------------------
   */
  cerr << " begin update... " << endl;
  TIME_STAGE(this->updateCells(), "Update cells", "4(b)");
  TIME_STAGE(this->updateECMManagers(), "Update ECM Managers", "4(c)");
  TIME_STAGE(this->updatePatches(), "Update Patches", "4(d)");
  TIME_STAGE(this->updateChem(), "Update chem", "4(a)");
#else // PROFILE_MAJOR_STEPS

  // For testing purposes:
  // cout << "  TNF: " << this->WorldChem.totalTNF << ", TGF: " <<
  // this->WorldChem.totalTGF << ", IL1beta: " << this->WorldChem.totalIL1beta
  // << endl;

  /* --------------------------- CHEMICAL DIFFUSION ---------------------------
   */
  this->diffuseCytokines();

  /* ------------------------------ CELL SEEDING ------------------------------
   */
  // this->seedCells(hours);

  /* ------------------------------ CELL FUNCTION -----------------------------
   */
  this->executeCells();

  /* ------------------------------ ECM FUNCTION ------------------------------
   */
  this->executeECMs();
  this->requestECMfragments();

#ifdef MODEL_SCAFFOLD
  /* ------------------------- UPDATE CaAlg Properties ------------------------
   */
  this->updateSwellingRatio();
  this->updateMassLoss();
  this->updateE();
#endif

  /* ----------------------- ATTRIBUTES SYNCHRONIZATION -----------------------
   */
  cerr << " begin update... " << endl;
  this->updateCells();
  this->updateECMManagers();
  this->updatePatches();
  this->updateChem();

#endif
  return 0;
}

/* -------------------------------------------------------------------------- */
/*                      MAJOR SECTION SUBROUTINES - begin                     */
/* -------------------------------------------------------------------------- */
/*
 * Chemical diffusion tick - delegated to ChemicalEnvironment (see go()
 * ordering). PDE numerics run inside ChemicalEnvironment (Diffusion3D).
 */

ChemicalEnvironment *BMWorld::chemical_environment() {
  return chemical_environment_.get();
}

const ChemicalEnvironment *BMWorld::chemical_environment() const {
  return chemical_environment_.get();
}

float BMWorld::world_total_tnf() const {
  if (chemical_environment_)
    return chemical_environment_->total_tnf();
  return WorldChem.totalTNF;
}

float BMWorld::world_total_tgf() const {
  if (chemical_environment_)
    return chemical_environment_->total_tgf();
  return WorldChem.totalTGF;
}

float BMWorld::world_total_il1beta() const {
  if (chemical_environment_)
    return chemical_environment_->total_il1beta();
  return WorldChem.totalIL1beta;
}

double BMWorld::tick_interval_minutes() const {
  if (chemical_environment_)
    return chemical_environment_->tick_interval_minutes();
  return kTickIntervalMinutes;
}

float BMWorld::chem_concentration(SpeciesId species, int patch_index) const {
  if (!chemical_environment_)
    return 0.f;
  return chemical_environment_->concentration_at(patch_index, species);
}

void BMWorld::chem_add_secretion(SpeciesId species, int patch_index,
                                 float delta) const {
  if (chemical_environment_)
    chemical_environment_->accumulate_secretion(patch_index, species, delta);
}

float BMWorld::chem_concentration_channel(int concentration_channel,
                                          int patch_index) const {
  if (!chemical_environment_)
    return 0.f;
  return chemical_environment_->concentration_at_channel(patch_index,
                                                         concentration_channel);
}

float BMWorld::chemotaxis_at(int patch_index) const {
  if (!chemical_environment_)
    return 0.f;
  return chemical_environment_->chemotaxis_at(patch_index);
}

void BMWorld::diffuseCytokines() {
  if (chemical_environment_)
    chemical_environment_->run_diffusion_phase(tick_interval_minutes());
}

void BMWorld::runCells() {
  int cellsSize =
      cells.size(); /* This is only an upper bound on cell list size. It is NOT
                       an actual count of cells (some entries are NULL) */
#pragma omp parallel for
  for (int i = 0; i < cellsSize; i++) {
    Cell *cell = cells.getDataAt(i);
    if (!cell)
      continue;
    cell->cellFunction();
  }
}

void BMWorld::executeCells() {
#ifdef PROFILE_CELL_FUNC
  TIME_STAGE(this->runCells(), "Cell Function: Chondrocytes", "	");
#else
  cout << " execute cells " << endl;
  this->runCells();
#endif
}

void BMWorld::executeECMs() {
  cerr << " ECM function  " << endl;
  int numPatches = (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny;
#pragma omp parallel for
  for (int in = 0; in < numPatches; in++) {
    if (worldECM[in].empty[read_t] == false)
      this->worldECM[in].ECMFunction();
  }
}

void BMWorld::requestECMfragments() {
  if (BMWorld::highTNFdamage == true) {
    cout << " high TNF damage " << endl;
    BMWorld::highTNFdamage = false;
    for (int in = 0; in < (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny; in++) {
#ifndef CALIBRATION
      if (this->chem_concentration(TNF, in) > 10.f) {
#else
      if (this->chem_concentration(TNF, in) > 10.f) {
#endif
        cout << " Degrade ECM " << endl;
        this->worldECM[in].fragmentNCollagen();
        this->worldECM[in].fragmentNAggrecan();
      }
    }
  }
}

void BMWorld::updateChemCPU() {
  if (!chemical_environment_)
    return;
  chemical_environment_->merge_and_reset_secretion();
  chemical_environment_->copy_totals_to(this->WorldChem);
}

void BMWorld::updateChem() { updateChemCPU(); }

void BMWorld::executeAllECMUpdates() {
  for (int iz = 0; iz < nz; iz++) {
#pragma omp parallel for
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        this->worldECM[in].updateECM();
      }
    }
  }
}

void BMWorld::executeAllECMResetRequests() {
  for (int iz = 0; iz < nz; iz++) {
#pragma omp parallel for
    for (int iy = 0; iy < ny; iy++) {
      for (int ix = 0; ix < nx; ix++) {
        int in = ix + iy * nx + iz * nx * ny;
        this->worldECM[in].resetrequests();
      }
    }
  }
}

void BMWorld::updateECMManagers() {
#ifdef PROFILE_ECM_UPDATE
  TIME_STAGE(this->executeAllECMUpdates(), "	updateECM()", "	");
  TIME_STAGE(this->executeAllECMResetRequests(), "	resetrequests()", "	");
#else
  this->executeAllECMUpdates();
  this->executeAllECMResetRequests();
#endif
}

/*
 * Steps:
 * 1. Perform updates
 * 2. Remove all dead cells
 * 3. If OMP, add cells from thread-local lists to corresponding global lists
 */
void BMWorld::updateCells() {
  cout << " previous tick's cells " << prevCells << endl;
  cerr << "	removing dead cells" << endl;
  int cellsSize = cells.size();
  liveCells = 0;
  int dcells = 0;
  deletedCells = 0;
#pragma omp parallel for
  for (int i = 0; i < cellsSize; i++) {
#ifdef _OMP
    int tid = omp_get_thread_num();
#else
    int tid = DEFAULT_TID;
#endif
    // Get pointer of cell i from the array chain
    Cell *cell = cells.getDataAt(i);
    if (!cell)
      continue; // cell was deleted
    cell->updateAgent();

    if (cell->isRealDead() == true) {
      dcells++;
    }

    if (cell->isAlive() == true) {
      liveCells++;
      // Update cell stage counts
      /* Added by MM to check types of cell stages and add to respective
       * counters: */
      if (typeid(*cell) == typeid(Stem)) {
        Stem::numOfStem++;
      } else if (typeid(*cell) == typeid(Progen)) {
        Progen::numOfProgen++;
      } else if (typeid(*cell) == typeid(NP)) {
        NP::numOfNP++;
      }

      // Remove dead cells
    } else if (cell->isAlive() == false) {
      // Get residing patch index and update its occupancy
      int in = cell->getIndex();
      this->worldPatch[in].clearOccupied();
      this->worldPatch[in].occupiedby[write_t] = nothing;
      this->worldPatch[in].dirty = true;
      /* Added by MM to check types of cell stages and subtract from respective
       * counters: */
      if (typeid(*cell) == typeid(Stem)) {
        Stem::numOfStem--;
      } else if (typeid(*cell) == typeid(Progen)) {
        Progen::numOfProgen--;
      } else if (typeid(*cell) == typeid(NP)) {
        NP::numOfNP--;
      }
      cells.deleteData(i, tid);
      delete cell;
      deletedCells++;
    }
  }
  // if (BMWorld::clock == 0) {
  //	prevCells = this->initialCells[0];
  // }
  // else if (BMWorld::clock > 0) {
  //	prevCells = liveCells;
  // }

  // deadCells += dcells;
  if (prevCells - cells.actualSize() >= 0) {
    deadCells += prevCells - cells.actualSize();
  } else if (prevCells - cells.actualSize() < 0) {
    deadCells += 0;
  }
  // deadCells += prevCells - cells.actualSize();
  Cell::numOfCells = cells.actualSize();
  // cout << " number of dead cells in this tick " << prevCells -
  // cells.actualSize() << endl;
  cout << " number of dead cells in this tick (prev - cells actual size) "
       << prevCells - cells.actualSize() << endl;
  cout << " number of dead cells in this tick (dcells) " << dcells << endl;
  cout << " number of cells now (actualSize) = " << cells.actualSize() << endl;
  cout << " number of live cells = " << liveCells << endl;
  cout << " total number of dead cells " << deadCells << endl;
  cout << " cell viability " << std::setprecision(3)
       << (liveCells / (liveCells + deadCells)) * 100 << "%" << endl;

// Add new cells
#ifdef _OMP
  /* In OMP, cells were only added to each thread's local list when
   * sproutAgent() was called. Thus, this step is needed to add those cells onto
   * the global lists. */
  cerr << "	updateCell() _OMP" << endl;
  // TODO: parallelize
  // int numThreads = omp_get_num_threads();
  int numThreads = std::max(atoi(std::getenv("OMP_NUM_THREADS")), 1);
  // cout << "		numThreads = " << numThreads << endl;
  for (int tid = 0; tid < numThreads; tid++) {
    // Cells
    vector<Cell *> *fvec_ptr = localNewCells[tid];
    for (vector<Cell *>::iterator cell_it = fvec_ptr->begin();
         cell_it != fvec_ptr->end(); cell_it++) {
      Cell *newCell = *cell_it;
      if (!cells.addData(newCell, tid)) {
        cerr << "Error: Could not add cell" << endl;
        exit(-1);
      }
    }
    fvec_ptr->clear();
  }
#endif
  prevCells = cells.actualSize();
}

/****************************************************************
 * MAJOR SECTION SUBROUTINES - end                              *
 ****************************************************************/
// NOTE: only use this function to sprout new_coll, new_agg in initialization.

void BMWorld::sproutAgent(int num, int patchType, int agentType, int xmin,
                          int xmax, int ymin, int ymax, int zmin, int zmax) {
#ifdef OPT_CELL_SEEDING
  if (xmin != 0 || xmax != nx || ymin != 0 || ymax != ny || zmin != 0 ||
      zmax != nz)
    sproutAgentInArea(num, patchType, agentType, xmin, xmax, ymin, ymax, zmin,
                      zmax);
  else if (patchType == CaAlg)
    sproutAgentInArea(num, patchType, agentType, xmin, xmax, ymin, ymax, zmin,
                      zmax);
  else
    sproutAgentInWorld(num, patchType, agentType);
#else
  // Target a specific area of the world
  sproutAgentInArea(num, patchType, agentType, xmin, xmax, ymin, ymax, zmin,
                    zmax);
#endif
}

void BMWorld::sproutAgentInArea(int num, int patchType, int agentType, int xmin,
                                int xmax, int ymin, int ymax, int zmin,
                                int zmax) {
  int count = 0;
  vector<int> patchlist;
  int *reservoir = new int[num];
  for (int i = 0; i < num; i++)
    reservoir[i] = -1;
  Patch *tempPatchPtr;
  Agent *tempAgentPtr;
  int in, agentIndex, max;
  for (int izz = zmin; izz < zmax + 1; izz++) {
    for (int iyy = ymin; iyy < ymax + 1; iyy++) {
      for (int ixx = xmin; ixx < xmax + 1; ixx++) {
        in = ixx + iyy * nx + izz * nx * ny;

        // Try another patch if this one is out of bounds or the wrong type or
        // occupied
        if (ixx < 0 || ixx >= nx || iyy < 0 || iyy >= ny || izz < 0 ||
            izz >= nz)
          continue;
        if (BMWorld::worldPatch[in].type[read_t] != patchType)
          continue;
        if (this->worldPatch[in].isOccupied() == false)
          patchlist.push_back(in);
      }
    }
  }
  for (int i = 0; i < num; i++) {
    if (patchlist.size() == 0) { // No available patches
      cout << " sprout agent error, no available patch within bounds! " << endl;
      delete[] reservoir;
      return;
    }
    int randnumber = rand() % patchlist.size();
    reservoir[i] = patchlist[randnumber]; // Prepare 'num' random patches
  }

  // Sprout agent on each patch in reservoir
  for (int i = 0; i < num; i++) {
    int in = reservoir[i];
    if (in < 0 || in > (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny)
      continue;
    switch (agentType) {
    case stem: {
      // cout << "new stem added" << endl; //added for debugging
      tempPatchPtr = &(this->worldPatch[in]);
      // std::shared_ptr<Cell> cell = std::make_shared<Stem>(tempPatchPtr);
      Stem *newStem = new Stem(tempPatchPtr);
#ifdef _OMP
      int tid = omp_get_thread_num();
      this->localNewCells[tid]->push_back(newStem);
#else
      if (!this->cells.addData(newStem, DEFAULT_TID)) {
        cerr << "Error: Could not add stem cell in sproutAgent()" << endl;
        exit(-1);
      }
#endif
      ///* Added by MM to check types of cell stages and add to respective
      /// counters: */
      // if (typeid(*this) == typeid(Stem)) {
      //	Stem::numOfStem++;
      // }
      // Cell::numOfCells++;

      this->worldPatch[in].setOccupied();
      this->worldPatch[in].occupiedby[write_t] = stem;
      this->worldPatch[in].dirty = true;
      // cout << "patch index " << in << endl; //added for debugging

      break;
    }
    case progen: {
      // cout << "new progen added" << endl; //added for debugging
      tempPatchPtr = &(this->worldPatch[in]);
      // std::shared_ptr<Cell> cell = std::make_shared<Progen>(tempPatchPtr);
      Progen *newProgen = new Progen(tempPatchPtr);
#ifdef _OMP
      int tid = omp_get_thread_num();
      this->localNewCells[tid]->push_back(newProgen);
#else
      if (!this->cells.addData(newProgen, DEFAULT_TID)) {
        cerr << "Error: Could not add pre-np cell in sproutAgent()" << endl;
        exit(-1);
      }
#endif
      ///* Added by MM to check types of cell stages and add to respective
      /// counters: */
      // if (typeid(*this) == typeid(Progen)) {
      //	Progen::numOfProgen++;
      // }
      // Cell::numOfCells++;

      this->worldPatch[in].setOccupied();
      this->worldPatch[in].occupiedby[write_t] = progen;
      this->worldPatch[in].dirty = true;

      break;
    }
    case np: {
      // cout << "new np added" << endl; //added for debugging
      tempPatchPtr = &(this->worldPatch[in]);
      // std::shared_ptr<Cell> cell = std::make_shared<NP>(tempPatchPtr);
      NP *newNP = new NP(tempPatchPtr);
#ifdef _OMP
      int tid = omp_get_thread_num();
      this->localNewCells[tid]->push_back(newNP);
#else
      if (!this->cells.addData(newNP, DEFAULT_TID)) {
        cerr << "Error: Could not add np cell in sproutAgent()" << endl;
        exit(-1);
      }
#endif
      ///* Added by MM to check types of cell stages and add to respective
      /// counters: */
      // if (typeid(*this) == typeid(NP)) {
      //	NP::numOfNP++;
      // }
      // Cell::numOfCells++;

      this->worldPatch[in].setOccupied();
      this->worldPatch[in].occupiedby[write_t] = np;
      this->worldPatch[in].dirty = true;

      break;
    }
    case orig_coll: {
      this->worldECM[in].ocollagen[write_t] =
          this->worldECM[in].ocollagen[read_t] + 1;
      this->worldECM[in].dirty = true;
      this->worldECM[in].isEmpty();
      break;
    }
    case orig_agg: {
      this->worldECM[in].oaggrecan[write_t] =
          this->worldECM[in].oaggrecan[read_t] + 1;
      this->worldECM[in].dirty = true;
      this->worldECM[in].isEmpty();
      break;
    }
    case new_coll: {
      this->worldPatch[in].initcollagen = true;
      this->worldECM[in].ncollagen[write_t] =
          this->worldECM[in].ncollagen[read_t] + 1;
      this->worldECM[in].dirty = true;
      this->worldECM[in].isEmpty();
      break;
    }
    case new_agg: {
      this->worldPatch[in].initaggrecan = true;
      this->worldECM[in].naggrecan[write_t] =
          this->worldECM[in].naggrecan[read_t] + 1;
      this->worldECM[in].dirty = true;
      this->worldECM[in].isEmpty();
      break;
    }
    }
  }
  delete[] reservoir;
  return;
}

#ifdef OPT_CELL_SEEDING
/*
 * Optimized by:
 *  - If sprout in tissue:
 *     (*) Randomly choosing a target patch:
 *          - If is tissue and occupied, sprout
 * 			    - Else repeat (*)
 *  - Else (sprout in blood):
 *     (**)	Look at the list of capillary patches initialized in the setup
 * stage
 *          - Create a list of unoccupied capillary patches
 *          - Pick randomly and sprout
 *          - Repeat until 'num' cells are sprouted
 */
void BMWorld::sproutAgentInWorld(int num, int patchType,
                                 int agentType) { // bool bloodORtiss
  int count = 0;
  vector<int> patchlist;
  int *reservoir = new int[num];
  for (int i = 0; i < num; i++)
    reservoir[i] = -1;
  Patch *tempPatchPtr;
  Agent *tempAgentPtr;
  int in, agentIndex, max;
  int totalNumPatches =
      this->(n - 1) x * this->(ny - 1) * this->(nz - 1); // Is this accurate?
  int numfound;
  int counter;
  int threshold;
  vector<int> unoccupiedCaps;
  switch (patchType) {
    for (int i = 0; i < num; i++) {
      if (reservoir[i] < 0 ||
          reservoir[i] > (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny)
        continue;
      switch (agentType) {
      case stem: {
        tempPatchPtr = &(this->worldPatch[reservoir[i]]);
        Stem *newStem = new Stem(tempPatchPtr);
#ifdef _OMP
        int tid = omp_get_thread_num();
        this->localNewCells[tid]->push_back(newStem);
#else
        if (!this->cells.addData(newStem, DEFAULT_TID)) {
          cerr << "Error: Could not add stem cell in sproutAgentInWorld()"
               << endl;
          exit(-1);
        }
#endif
        // this->worldPatch[reservoir[i]].occupied = true;
        this->worldPatch[reservoir[i]].setOccupied();
        this->worldPatch[reservoir[i]].occupiedby = stem;
        break;
      }
      case progen: {
        tempPatchPtr = &(this->worldPatch[reservoir[i]]);
        Progen *newProgen = new Progen(tempPatchPtr);
#ifdef _OMP
        int tid = omp_get_thread_num();
        this->localNewCells[tid]->push_back(newProgen);
#else
        if (!this->cells.addData(newProgen, DEFAULT_TID)) {
          cerr << "Error: Could not add pre-np cell in sproutAgentInWorld()"
               << endl;
          exit(-1);
        }
#endif
        // this->worldPatch[reservoir[i]].occupied = true;
        this->worldPatch[reservoir[i]].setOccupied();
        this->worldPatch[reservoir[i]].occupiedby = progen;
        break;
      }
      case np: {
        tempPatchPtr = &(this->worldPatch[reservoir[i]]);
        NP *newNP = new NP(tempPatchPtr);
#ifdef _OMP
        int tid = omp_get_thread_num();
        this->localNewCells[tid]->push_back(newNP);
#else
        if (!this->cells.addData(newNP, DEFAULT_TID)) {
          cerr << "Error: Could not add np cell in sproutAgentInWorld()"
               << endl;
          exit(-1);
        }
#endif
        // this->worldPatch[reservoir[i]].occupied = true;
        this->worldPatch[reservoir[i]].setOccupied();
        this->worldPatch[reservoir[i]].occupiedby = np;
        break;
      }
      }
    }
    delete[] reservoir;
    return;
  }
#endif // OPT_CELL_SEEDING

  int BMWorld::countPatchType(int whichType) {
    if (whichType == CaAlg) {
      Patch::numOfEachTypes[whichType] = 0;
      for (int iz = 0; iz < this->nz; iz++) {
        int currCount = 0;
#pragma omp parallel for reduction(+ : currCount)

        for (int iy = 0; iy < this->ny; iy++) {
          for (int ix = 0; ix < this->nx; ix++) {
            int in = ix + iy * nx + iz * nx * ny;
            if (this->worldPatch[in].type[read_t] == whichType) {
              currCount++;
              // Patch::numOfEachTypes [whichType]++;
            }
          }
        }
        Patch::numOfEachTypes[whichType] += currCount;
      }
    } else if (whichType == damage) {
      Patch::numOfEachTypes[whichType] = 0;
      for (int iz = 0; iz < this->nz; iz++) {
        int currCount = 0;
#pragma omp parallel for reduction(+ : currCount)
        for (int iy = 0; iy < this->ny; iy++) {
          for (int ix = 0; ix < this->nx; ix++) {
            int in = ix + iy * nx + iz * nx * ny;
            currCount += this->worldPatch[in].damage[read_t];
          }
        }
        Patch::numOfEachTypes[whichType] += currCount;
      }
    } else
      cout << "type must be 0, 1, 2 , 3 or 4!"
           << endl; // cout << "type must be 0, 1, 2 , 3, 4 or 5!" << endl;
    return Patch::numOfEachTypes[whichType];
  }

  int BMWorld::mmToPatch(double mm) { return mm * (this->patchpermm); }

  int BMWorld::reportTick(int hour, int day) { return (hour * 2 + day * 48); }

  double BMWorld::reportMinute() { return (BMWorld::clock) * 30; }

  double BMWorld::reportHour() { return (BMWorld::clock) / 2; }

  double BMWorld::reportDay() { return (BMWorld::clock) / 48; }

  int BMWorld::countNeighborPatchType(int ix, int iy, int iz, int patchType) {
    int neighborcount = 0;
    for (int dx = -1; dx < 2; dx++) {
      for (int dy = -1; dy < 2; dy++) {
        for (int dz = -1; dz < 2; dz++) {
          int neighborindex = (ix + dx) + (iy + dy) * nx + (iz + dz) * ny * nx;
          if (ix + dx < 0 || ix + dx >= nx || iy + dy < 0 || iy + dy >= ny ||
              iz + dz < 0 || iz + dz >= nz)
            continue;
          if (dx == 0 && dy == 0 && dz == 0)
            continue;
          if (Agent::agentPatchPtr[neighborindex].type[read_t] == patchType)
            neighborcount++;
        }
      }
    }
    return neighborcount;
  }

#ifdef MODEL_SCAFFOLD
  void BMWorld::updateSwellingRatio() {
    float Alg_ww = this->Alg_wv / (this->Alg_wv);
    double tmin = reportMinute();

/* Calculate Swelling Ratio (Q) given Alg concentration (% w/w) at a given time
 * (in minutes) Q = (a * Alg_ww +  b) * t_min + (c * Alg_ww + d)
 *
 *       Swelling ratio favorable for cell adhesion, growth, diffusion of
 * nutrients Depends on Alg content (% w/w) Rapid swelling in initial 10 min,
 * with slight increase until ~24h
 */
#ifdef CALIBRATION
    // this->Q = (this->SwellRatio[0]*Alg_ww + this->SwellRatio[1])*log(tmin) +
    // (this->SwellRatio[2]*Alg_ww + this->SwellRatio[3]); old
    BMWorld::SwellRatio[0] - BMWorld::SwellRatio[1] * (this->reportDay()) -
        BMWorld::SwellRatio[2] * (this->Alg_wv) -
        BMWorld::SwellRatio[3] * (this->reportDay()) * (this->pXL) +
        BMWorld::SwellRatio[4] * (this->Alg_wv) * (this->pXL);
#else
    this->Q = (0.4 * Alg_ww + 0.4) * log(tmin) + (3 * Alg_ww + 7.9);
#endif
    /* Swelling reduces effective diffusivity: D_eff = base_D * Q in
     * SpeciesRegistry. */
    if (this->Q > 0.0f && chemical_environment_)
      chemical_environment_->set_swelling_ratio(static_cast<double>(this->Q));
    // this->Q=0;
    // cout << " this->Q =" << (this->SwellRatio[0]<<"*"<<Alg_ww<< " + "<<
    // this->SwellRatio[1])<<"*log("<<tmin<<") +
    // ("<<this->SwellRatio[2]<<"*"<<Alg_ww<< " + "<<this->SwellRatio[3]<<")" <<
    // endl; cout << " Swelling Ratio: " << this->Q << endl;
  }

  void BMWorld::updateMassLoss() {
    float Alg_ww = this->Alg_wv / (this->Alg_wv);
    float tweek = reportDay() / 7;
    float w_t;

/* Calculate Mass Loss (w_t) of gel with given Alg concentration (% w/w) at
 * current time (in weeks) w_t = (a * Alg_ww +  b) * t_weeks + (c * Alg_ww + d)
 *
 *       Degree of dregradation depends greatly on Alg (% w/w) content
 *       Highest weight loss percentage occurs during first week in vitro, with
 * little weight loss over next 3 weeks
 */
#ifdef CALIBRATION
    w_t = this->MassLoss[0] + this->MassLoss[1] * (pXL) +
          this->MassLoss[2] * (reportDay()) -
          this->MassLoss[3] * (pXL) * (reportDay());
#else
    w_t =
        0.234 + 7.785 * (pXL) + 0.15 * (reportDay()) -
        1.36 * (pXL) *
            (reportDay()); //(17.6*Alg_ww - 0.9)*log(tweek) + (60*Alg_ww + 5.3);
#endif
    if (w_t < 0)
    w_t = 0; // no negative mass loss

    // If there is % mass loss since last call, "degrade" % CaAlg patches and
    // replace with tissue
    if (w_t > this->w) {
      float changeInPatches = (0.01) * (w_t - this->w) * BMWorld::initialCaAlg;
      this->degradeCaAlg(changeInPatches);
    }

    this->w = w_t;
    // cout << " w_t = "<<this->MassLoss[0]<<" + "<<
    // this->MassLoss[1]<<"*"<<(pXL)<<" +
    // "<<this->MassLoss[2]<<"*"<<(reportDay()) <<" - "
    // <<this->MassLoss[3]<<"*"<<(pXL)<<"*"<<(reportDay()) << endl; cout << "
    // Mass Loss (%): " << this->w << endl; cout << " Number of Ca-Alg patches:
    // " << this->countPatchType(CaAlg) << endl;
  }
  
  void BMWorld::updateE() {
	  BMWorld::E = BMWorld::E_inf + (BMWorld::E_0 - BMWorld::E_inf) * exp(-(BMWorld::clock * 30 * 60) / BMWorld::t); // converts tick to seconds
  }
#endif // MODEL_SCAFFOLD

#ifdef MODEL_SCAFFOLD
  void BMWorld::degradeCaAlg(int numOfPatches) {
    int xmin = 0;
    int xmax = nx;
    int ymin = 0;
    int ymax = ny;
    int zmin = 0;
    int zmax = nz;
    int patchType = CaAlg;
    vector<int> patchlist;
    int *reservoir = new int[numOfPatches];
    for (int i = 0; i < numOfPatches; i++)
      reservoir[i] = -1;
    int in, agentIndex, max;
    int count = 0;

    // Make list of possible CaAlg Patches to degrade
    for (int iz = zmin; iz < zmax; iz++) {
      for (int iy = ymin; iy < ymax; iy++) {
        for (int ix = xmin; ix < xmax; ix++) {
          in = ix + iy * nx + iz * nx * ny;

          // Try another patch if this one is out of bounds or the wrong type or
          // occupied
          if (ix < 0 || ix >= nx || iy < 0 || iy >= ny || iz < 0 || iz >= nz)
            continue;
          if (BMWorld::worldPatch[in].type[read_t] != patchType)
            continue;
          patchlist.push_back(in);
        }
      }
    }

    // Choose random patches from patch list
    for (int i = 0; i < numOfPatches; i++) {
      if (patchlist.size() == 0) { // No available patches
        // cout << " CaAlg degrade error, no available patch within bounds! " <<
        // endl;
        delete[] reservoir;
        return;
      }
      int randnumber = rand() % patchlist.size();
      reservoir[i] = patchlist[randnumber]; // Prepare 'num' random patches
    }

    // Degrade 'numOfPatches" number of CaAlg patches in reservoir list
    for (int i = 0; i < numOfPatches; i++) {
      int in = reservoir[i];
      if (in < 0 || in > (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny)
        continue;
      BMWorld::worldPatch[in].type[write_t] = nothing;
      BMWorld::worldPatch[in].color[write_t] = cnothing;
      BMWorld::worldPatch[in].dirty = true;
      count++;
    }
    delete[] reservoir;
  }
#endif // MODEL_SCAFFOLD

  void BMWorld::debugInfo() {
    int alive = 0;
    int dead = 0;
    int stemSize = 0;
    int progenSize = 0;
    int npSize = 0;
    int cellsSize = cells.size();
    for (int i = 0; i < cellsSize; i++) {
      Cell *cell = cells.getDataAt(i);
      if (!cell)
        continue;
      if (cell->isAlive() == false) {
        dead++;
      } else if (cell->isAlive() == true) {
        alive++;
      }
      // if (cell->activate[read_t] == false) f++;
      else if (typeid(*cell) == typeid(Stem)) {
        stemSize++;
      } else if (typeid(*cell) == typeid(Progen)) {
        progenSize++;
      } else if (typeid(*cell) == typeid(NP)) {
        npSize++;
      }
      // else af++;
    }

    int numCaAlg = 0;
    numCaAlg = countPatchType(CaAlg);
    cout << " total patches: " << numCaAlg << endl;
    cout << " alive cells: " << alive << endl;
    cout << " dead cells: " << dead << endl;
    cout << " total cells: " << cells.actualSize() << endl;
    cout << " stem cells: " << stemSize << endl;
    cout << " pre-np cells: " << progenSize << endl;
    cout << " np cells: " << npSize << endl;
  }

  /*
   * Steps:
   *   1. If OMP, add cells from thread-local lists to corresponding global
   * lists
   *   2. Perform updates
   */
  void BMWorld::updateCellsInitial() {
// Cell lists should be empty
// Add new cells
#ifdef _OMP
    cerr << "	updateCell() _OMP" << endl;
    // TODO: parallelize

    int numThreads = omp_get_num_threads();
    for (int tid = 0; tid < numThreads; tid++) {
      /* ------------------------------ Cells ------------------------------ */
      vector<Cell *> *fvec_ptr = localNewCells[tid];
      for (vector<Cell *>::iterator cell_it = fvec_ptr->begin();
           cell_it != fvec_ptr->end(); cell_it++) {
        Cell *newCell = *cell_it;
        if (!cells.addData(newCell, tid)) {
          cerr << "Error: Could not add cell" << endl;
          exit(-1);
        }
      }
      fvec_ptr->clear();
    }
#endif

    /* ------------------------------ Cells ------------------------------ */
    // No need for deletion since these are new cells
    int cellsSize = cells.size();
#pragma omp parallel for
    for (int i = 0; i < cellsSize; i++) {
#ifdef _OMP
      int tid = omp_get_thread_num();
#else
    int tid = DEFAULT_TID;
#endif
      Cell *cell = cells.getDataAt(i);
      if (!cell)
        continue;
      cell->updateAgent();
      /* Added by MM to check types of cell stages and add to respective
       * counters: */
      if (typeid(*cell) == typeid(Stem)) {
        Stem::numOfStem++;
      } else if (typeid(*cell) == typeid(Progen)) {
        Progen::numOfProgen++;
      } else if (typeid(*cell) == typeid(NP)) {
        NP::numOfNP++;
      }
    }
    Cell::numOfCells = cells.actualSize();
  }

  int BMWorld::userInput() {
    // Read input parameters from user-specified file
    ifstream infile(util::getInputFileName());

    int numChem = -1;
    // TODO: Make this check for specific tag (field name)

    if (infile.is_open()) {
      char garbage[100];
      /* -------------------------------- CHEMICALS
       * ------------------------------- */
      // cout << "Reading the number of baseline chemicals..." << endl;
      float temp;
      infile >> garbage;
      infile >> temp;
      // cout << garbage;
      // cin >> garbage;

      numChem = temp;
      this->baselineChem.resize(temp);
      cout << "The number of baseline chemicals are: "
           << this->baselineChem.size() << endl;
      for (int ichem = 0; ichem < baselineChem.size(); ichem++) {
        // cout << "Reading the baselineChemical " << ichem << endl;
        infile >> garbage;
        infile >> this->baselineChem[ichem];
        // cout << garbage;
        // cin >> garbage;
        cout << "baselineChem " << ichem << " is " << this->baselineChem[ichem]
             << endl;
      }

      /* ---------------------------------- CELLS
       * --------------------------------- */
      // cout << "Reading the initial number of types of cells..." << endl;
      int tempCells;
      infile >> garbage;
      infile >> tempCells;
      // cout << garbage;
      // cin >> garbage;
      this->initialCells.resize(tempCells);
      cout << "The number of types of cells is: " << this->initialCells.size()
           << endl;
      for (int icell = 0; icell < initialCells.size(); icell++) {
        // cout << "Reading the initial # of cell type " << icell + 1 << endl;
        infile >> garbage;
        infile >> this->initialCells[icell];
        // cout << garbage;
        // cin >> garbage;
        cout << "# of cell " << icell << " is " << this->initialCells[icell]
             << endl;
      }

      /* -------------------------- Ca-Alg PROPERTIES --------------------------
       */
      /* NOTE FOR STEM CELL IVDBM-ABM: none of this is being used anywhere */
      /* currently. Commented out and replaced with the important params that go
       */
      /* into the elasticity equation: alginate kDas and calcium concentration
       */
      /* in mM. */

      ////cout << "Reading 1% (w/v) Alg volume (mL)" << endl;
      // infile >> garbage;
      // infile >> this->Alg_v;
      // cout << "Volume of Alg (mL) = " << this->Alg_v << endl;

      ////cout << "Reading 1.67% (w/v) Ca volume (mL)" << endl;
      // infile >> garbage;
      // infile >> this->Ca_v;
      // cout << "Volume of Ca (mL) = " << this->Ca_v << endl;
      //
      // float totalVolume = this->Alg_v + this->Ca_v;

      // this->Alg_wv = 2;
      // this->Alg_wv = 1.95;

      // cout << "Total Volume from file (mL): " << totalVolume << endl;

      infile >> garbage;
      infile >> this->Alg_wv;
      cout << "Concentration of Alg w/v (%) = " << this->Alg_wv << endl;

      infile >> garbage;
      infile >> this->highMW_alg;
      cout << "Ratio component of high MW alginate (integer) = "
           << this->highMW_alg << endl;

      infile >> garbage;
      infile >> this->lowMW_alg;
      cout << "Ratio component of low MW alginate (integer) = "
           << this->lowMW_alg << endl;

      infile >> garbage;
      infile >> this->pXL;
      cout << "Concentration of Ca crosslinker (mM) = " << this->pXL << endl;
      
      infile >> garbage;
		  infile >> this->peptide;
		  cout << "Type of peptide conjugation = " << this->peptide << endl;

      /* --------------------------- CYTOKINE PROPERTIES
       * -------------------------- */
      /* Species diffusivity and baselines come from chemical_environment.json.
       */
      infile.close();
      cout << "-------------------------------------------" << endl;
    } // end of if file opens properly
    else
      cout << "Cannot open file!" << endl;

    // this->baselineChem[3]; //added manually for testing
    // cout << "The number of baseline chemicals are: " <<
    // this->baselineChem.size() << endl;

    // this->initialCells[3]; //added manually for testing
    // cout << "The number of types of cells is: " << this->initialCells.size()
    // << endl;

    // this->Alg_v = 0.0275; //added manually for testing
    // this->Ca_v = 0.005;

    // float totalVolume = this->Alg_v + this->Ca_v;

    // cout << "Total Volume (mL): " << totalVolume << endl;

    return 0;
  }

  // void BMWorld::outputWorld_csv() {
  //	if (this->clock == 0) {
  //		remove("output/Output_Biomarkers.csv");
  //		remove("output/tgf_line.csv");
  //
  //		ofstream output_file("output/Output_Biomarkers.csv", ios::app);
  //		ofstream tgf_file("output/tgf_line.csv", ios::app);
  //
  //		output_file << "clock (30 min)" << "," << "Day" << "," << "Total
  // TNF (pg)" << "," << "Total IL1b (pg)" << "," << "Total TGF (pg)" << "," <<
  //"Collagen (ug)" << "," << "Aggrecan (ug)" << "," << "Total Cells" << "," <<
  //"Stem Cells" << ", Pre-NP Cells" << ", NP Cells" << ", Live Cells" << ",
  // Dead Cells" << ", Elastic Modulus(kPa) " << ", Swelling Ratio " << ", Mass
  // Loss(%) " << ", Alginate_wv(%)" << ", Alginate_Mw(kDa)" << ", Ca_XL(M)" <<
  //", Viability Rate(%)" << ", Differentiation (%)" << endl; //output_file <<
  //"Tropocollagen" << ", " << "Collagen" << ", " << "FragentedCollagen" << ", "
  //<< "Tropoaggrecan" << ", " << "Aggrecan" << ", " << "FragmentedAggrecan" <<
  //", " << "HA" << ", " << "FragmentedHA" << ", " << "Damage" endl;
  //		output_file.close();
  //
  //		// Header: patch x-indices
  //		for (int xi = 0; xi <= nx / 2; xi++) {
  //			//o2_file << "x=" << xi << (xi < nx / 2 ? "," : "\n");
  //			tgf_file << "x=" << xi << (xi < nx / 2 ? "," : "\n");
  //		}
  //		//o2_file.close();
  //		tgf_file.close();
  //	}
  //
  //	ofstream output_file("output/Output_Biomarkers.csv", ios::app);
  //	ofstream tgf_file("output/tgf_line.csv", ios::app);
  //
  //	int f = 0; int af = 0;
  //	int orig_coll = 0; int frag_coll = 0; double new_coll = 0;
  //	int orig_agg = 0; int frag_agg = 0; double new_agg = 0;
  //	int HA = 0; int fHA = 0;
  //	int stemSize = 0; int progenSize = 0; int npSize = 0;
  //	int alive = 0; int dead = 0;
  //	float cellViability = 0;
  //	float perDiff = 0; // percent of NP cells out of total cells at a given
  // tick 	int x = 0;
  //
  //	int cellsSize = cells.size();
  //	for (int i = 0; i < cellsSize; i++) {
  //		Cell* cell = cells.getDataAt(i);
  //		if (!cell) continue;
  //		if (cell->isAlive() == false) continue;
  //		if (cell->activate[read_t] == false) f++;
  //		if (typeid(*cell) == typeid(Stem)) {
  //			stemSize++;
  //		}
  //		else if (typeid(*cell) == typeid(Progen)) {
  //			progenSize++;
  //		}
  //		else if (typeid(*cell) == typeid(NP)) {
  //			npSize++;
  //		}
  //		else af++;
  //	}
  //
  //	cout << " total cells: " << cells.actualSize() << endl;
  //	cout << " stem cells: " << stemSize << endl;
  //	cout << " pre-np cells: " << progenSize << endl;
  //	cout << " np cells: " << npSize << endl;
  //
  //	for (int in = 0; in < (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny;
  // in++) { 		orig_coll += this->worldECM[in].ocollagen[read_t];
  // new_coll += this->worldECM[in].ncollagen[read_t]; 		frag_coll +=
  // this->worldECM[in].fcollagen[read_t]; 		orig_agg +=
  // this->worldECM[in].oaggrecan[read_t]; 		new_agg +=
  // this->worldECM[in].naggrecan[read_t]; 		frag_agg +=
  // this->worldECM[in].faggrecan[read_t];
  //	}
  //
  //	//calculating viable cells
  //	//for (int i = 0; i < cellsSize; i++) {
  //	//	Cell* cell = cells.getDataAt(i);
  //	//	if (!cell) continue;
  //	//	if (cell->isAlive() == true) {
  //	//		alive++;
  //	//	}
  //	//	else if (cell->isAlive() == false) {
  //	//		dead++;
  //	//	}
  //	//}
  //	//cellViability = (static_cast<float>(liveCells) / (liveCells +
  //(prevCells - cells.actualSize()))) * 100; 	if (this->clock == 0) {
  //		cellViability = 100.0;
  //	}
  //	else {
  //		cellViability = (static_cast<float>(liveCells) / (liveCells +
  // deadCells)) * 100;
  //	}
  //	perDiff = (static_cast<float>(npSize) / cells.actualSize()) * 100;
  //
  //	this->countPatchType(damage);
  //	output_file << this->clock << ",";
  //	output_file << (this->clock) / 48 << ",";
  //	output_file << this->WorldChem.totalTNF << ",";
  //	output_file << this->WorldChem.totalIL1beta << ",";
  //	output_file << this->WorldChem.totalTGF << ",";
  //	output_file << fixed << std::setprecision(5) << new_coll << "," <<
  // new_agg << ","; //ECM 	//output_file << orig_coll << "," << new_coll <<
  //"," << frag_coll << "," << orig_agg << "," ; 	//output_file << new_agg
  //<< "," << frag_agg << "," << HA << "," << fHA << "," <<
  // Patch::numOfEachTypes[4] << "," ;
  //	//output_file << fixed << std::setprecision(5) << new_coll << "," <<
  // new_agg << ","; //ECM 	//output_file << orig_coll << "," << new_coll <<
  //"," << frag_coll << "," << orig_agg << "," ; 	//output_file << new_agg
  //<< "," << frag_agg << "," << HA << "," << fHA << "," <<
  // Patch::numOfEachTypes[4] << "," ;
  //
  //	//output_file << af << "," << f+af << ","; //cells
  //	output_file << cells.actualSize() << "," << stemSize << "," <<
  // progenSize << "," << npSize << "," << liveCells << "," << deadCells << ",";
  //// cell counts
  //
  // #ifdef MODEL_SCAFFOLD
  //	output_file << this->E << " , " << this->Q << ", " << this->w << "," <<
  // this->Alg_wv << "," << this->Alg_Mn << "," << this->pXL << "," <<
  // cellViability << "," << perDiff << endl; #else 	output_file << endl;
  // #endif 	output_file.close();
  //
  //	cout << " Collagen: " << new_coll << endl;
  //	cout << " Aggrecan: " << new_agg << endl;
  //
  //	tgf_file << fixed << std::setprecision(10);
  //
  //	// Record chem along the line - one row per tick
  //	for (int xi = 0; xi <= nx / 2; xi++) {
  //		//o2_file << o2Line[xi] << (xi < nx / 2 ? "," : "\n");
  //		tgf_file << tgfLine[xi] << (xi < nx / 2 ? "," : "\n");
  //	}
  //
  //	tgf_file.close();
  //
  //	prevCells = cells.actualSize();
  // }
  char* BMWorld::get_output_filename() {
    return util::outputFileName;
    //return "output/Output_Biomarkers.csv";
  }

  vector<string> BMWorld::get_agent_type_names() {
    return {"Stem", "Progen", "NP"};
  }

  void BMWorld::count_agent_types(map<string, int> & agent_counts) {
    int cellsSize = cells.size();
    for (int i = 0; i < cellsSize; i++) {
      Cell *cell = cells.getDataAt(i);
      if (!cell || !cell->isAlive())
        continue;
      if (typeid(*cell) == typeid(Stem))
        agent_counts["Stem"]++;
      else if (typeid(*cell) == typeid(Progen))
        agent_counts["Progen"]++;
      else if (typeid(*cell) == typeid(NP))
        agent_counts["NP"]++;
    }
  }

  int BMWorld::get_total_agent_count() { return cells.actualSize(); }

  vector<string> BMWorld::get_env_type_names() {
    return {"ncollagen", "naggrecan"};
  }

  void BMWorld::count_env(map<string, float> & env_counts) {
    for (int in = 0; in < (nx - 1) + (ny - 1) * nx + (nz - 1) * nx * ny; in++) {
      env_counts["ncollagen"] += this->worldECM[in].ncollagen[read_t];
      env_counts["naggrecan"] += this->worldECM[in].naggrecan[read_t];
    }
  }

  void BMWorld::write_data_row(std::ofstream & file,
                               std::map<std::string, int> & agent_counts,
                               std::map<std::string, float> & env_counts) {

    // biomaterial world-specific chemicals
    file << this->WorldChem.totalTNF << "," << this->WorldChem.totalIL1beta
         << "," << this->WorldChem.totalTGF << ",";

    // ecm types
    file << fixed << setprecision(5) << env_counts["ncollagen"] << ","
         << env_counts["naggrecan"] << ",";

    // agent counts
    file << get_total_agent_count() << "," << agent_counts["Stem"] << ","
         << agent_counts["Progen"] << "," << agent_counts["NP"] << ","
         << liveCells << "," << deadCells << ",";

    // viability and differentiation- specific to this biomaterial world
    float viability = calculate_viability();
    float perDiff = calculate_pct_differentiated(agent_counts);

    // scaffold-specific columns
#ifdef MODEL_SCAFFOLD
    file << this->E << "," << this->Q << "," << this->w << "," << this->Alg_wv
         << "," << this->Alg_Mn << "," << this->pXL << "," << viability << ","
         << perDiff << endl;
#else
  file << viability << "," << perDiff << endl;
#endif
  }

  // private helpers
  float BMWorld::calculate_viability() {
    if (this->clock == 0)
      return 100.0f;
    return (static_cast<float>(liveCells) / (liveCells + deadCells)) * 100;
  }

  float BMWorld::calculate_pct_differentiated(map<string, int> & agent_counts) {
    return (static_cast<float>(agent_counts["NP"]) / get_total_agent_count()) *
           100;
  }

  void BMWorld::write_csv_header(ofstream & file) {
    file << "clock (30 min)" << "," // written by skeleton
         << "Day" << ","            // written by skeleton
         << "Total TNF (pg)" << "," // everything below written by this hook
         << "Total IL1b (pg)" << ","
         << "Total TGF (pg)" << ","
         << "Collagen (ug)" << ","
         << "Aggrecan (ug)" << ","
         << "Total Cells" << ","
         << "Stem Cells" << ","
         << "Pre-NP Cells" << ","
         << "NP Cells" << ","
         << "Live Cells" << ","
         << "Dead Cells" << ","
         << "Elastic Modulus(kPa)" << ","
         << "Swelling Ratio" << ","
         << "Mass Loss(%)" << ","
         << "Alginate_wv(%)" << ","
         << "Alginate_Mw(kDa)" << ","
         << "Ca_XL(M)" << ","
         << "Viability Rate(%)" << ","
         << "Differentiation (%)" << endl;
  }

  // extra output - for IVDBM-ABM,
  // a TGF line for measuring TGF for
  // diffusion debugging purposes
  void BMWorld::write_auxiliary_header() {
    remove("output/tgf_line.csv");
    ofstream tgf_file("output/tgf_line.csv", ios::app);
    for (int xi = 0; xi <= nx / 2; xi++)
      tgf_file << "x=" << xi << (xi < nx / 2 ? "," : "\n");
    tgf_file.close();
  }

  void BMWorld::write_auxiliary_outputs() {
    ofstream tgf_file("output/tgf_line.csv", ios::app);
    tgf_file << fixed << setprecision(10);
    for (int xi = 0; xi <= nx / 2; xi++)
      tgf_file << tgfLine[xi] << (xi < nx / 2 ? "," : "\n");
    tgf_file.close();
  }

  void BMWorld::patchassign_csv() {
    int in = 0;
    int Number = 0;
    for (int iz = 0; iz < nz; iz++) {
      char patchassign[50] = "output/patchassign";
      char cells[50] = "output/cells_read";
      char cells_w[50] = "output/cells_write";
      char initcollagen[50] = "output/initcollagen";
      char initaggrecan[50] = "output/initaggrecan";
      char initHA[50] = "output/initHA";
      char damagezone[50] = "output/damagezone";
      char initialdamage[50] = "output/initialdamage";
      char extension[10] = ".csv";
      char tempNumber[20] = "";
      sprintf(tempNumber, "_t%3.0f_z%d", this->clock, Number);
      strcat(patchassign, tempNumber);
      strcat(patchassign, extension);
      strcat(cells, tempNumber);
      strcat(cells, extension);
      strcat(cells_w, tempNumber);
      strcat(cells_w, extension);
      strcat(initHA, tempNumber);
      strcat(initHA, extension);
      strcat(initcollagen, tempNumber);
      strcat(initcollagen, extension);
      strcat(initaggrecan, tempNumber);
      strcat(initaggrecan, extension);
      strcat(damagezone, tempNumber);
      strcat(damagezone, extension);
      strcat(initialdamage, tempNumber);
      strcat(initialdamage, extension);

      // Patch Assign
      ofstream output_file(patchassign, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].type[read_t] == damage ||
              worldPatch[in].damage[read_t] != 0) {
            output_file << "x";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file << "?";
          }
        }
        output_file << endl;
      }
      output_file.close();

      // initHA
      ofstream output_file1(initHA, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].initHA == true) {
            output_file1 << "u";
            continue;
          }
          if (worldPatch[in].type[read_t] == damage ||
              worldPatch[in].damage[read_t] != 0) {
            output_file1 << "x";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file1 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file1 << "?";
          }
        }
        output_file1 << endl;
      }
      output_file1.close();

      // Damage Zone
      ofstream output_file2(damagezone, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;

          if (worldPatch[in].inDamzone == true) {
            output_file2 << "z";
            continue;
          }
          if (worldPatch[in].type[read_t] == damage ||
              worldPatch[in].damage[read_t] != 0) {
            output_file2 << "x";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file2 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file2 << "?";
          }
        }
        output_file2 << endl;
      }
      output_file2.close();

      // Initial Damage
      ofstream output_file3(initialdamage, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldECM[in].oaggrecan[read_t] != 0 &&
              worldPatch[in].damage[read_t] != 0) {
            output_file3 << "g";
            continue;
          }
          if (worldECM[in].oaggrecan[read_t] != 0) {
            output_file3 << "m";
            continue;
          }
          if (worldPatch[in].type[read_t] == damage ||
              worldPatch[in].damage[read_t] != 0) {
            output_file3 << "x";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file3 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file3 << "?";
          }
        }
        output_file3 << endl;
      }
      output_file3.close();

      // initCollagen
      ofstream output_file4(initcollagen, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldECM[in].ocollagen[read_t] != 0) {
            output_file4 << "k";
            continue;
          }
          if (worldECM[in].fcollagen[read_t] != 0) {
            output_file4 << "f";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file4 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file4 << "?";
          }
        }
        output_file4 << endl;
      }
      output_file4.close();

      // initAggrecan
      ofstream output_file5(initaggrecan, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldECM[in].oaggrecan[read_t] != 0) {
            output_file5 << "m";
            continue;
          }
          if (worldECM[in].faggrecan[read_t] != 0) {
            output_file5 << "f";
            continue;
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file5 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file5 << "?";
          }
        }
        output_file5 << endl;
      }
      output_file5.close();

      // Cells
      ofstream output_file6(cells, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].isOccupied()) {
            if (worldPatch[in].occupiedby[read_t] == stem) {
              output_file6 << "f";
              continue;
            } else if (worldPatch[in].occupiedby[read_t] == progen) {
              output_file6 << "g";
              continue;
            } else if (worldPatch[in].occupiedby[read_t] == np) {
              output_file6 << "h";
              continue;
            }
          }
          if (worldPatch[in].type[read_t] == nothing) {
            output_file6 << "-";
          }
          if (worldPatch[in].type[read_t] == unidentifiable) {
            output_file6 << "?";
          }
        }
        output_file6 << endl;
      }
      output_file6.close();

      ofstream output_file7(cells_w, ios::app);
      for (int iy = 0; iy < ny; iy++) {
        for (int ix = 0; ix < nx; ix++) {
          in = ix + iy * nx + iz * nx * ny;
          if (worldPatch[in].isOccupiedWrite()) {
            if (worldPatch[in].occupiedby[read_t] == stem) {
              output_file6 << "f";
              continue;
            } else if (worldPatch[in].occupiedby[read_t] == progen) {
              output_file6 << "g";
              continue;
            } else if (worldPatch[in].occupiedby[read_t] == np) {
              output_file6 << "h";
              continue;
            }
          }
          if (worldPatch[in].type[write_t] == nothing) {
            output_file7 << "-";
          }
          if (worldPatch[in].type[write_t] == unidentifiable) {
            output_file7 << "?";
          }
        }
        output_file7 << endl;
      }
      output_file7.close();

      Number++;
    }
  }