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
#define	BMWORLD_H

#include "../World.h"
#include "../../FieldVariable/Usr_FieldVariables/Chemical.h"
#include "../../Agent/Usr_Agents/Cell.h"
#include "../../ECM/ECM.h"
#include "../../ArrayChain/ArrayChain.h"
#include "../../common.h"

#include <stdlib.h>
#include <vector>
#include <new>
#include <map>

#ifdef GPU_DIFFUSE
// Include CUDA runtime and CUFFT
#include <cuda_runtime.h>
#include <cufft.h>

// Helper functions for CUDA
//#include <helper_functions.h>
//#include <helper_cuda.h>

#include "../../Diffusion/convolutionFFT2D_common.h"
#endif	// GPU_DIFFUSE


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
 BMWORLD (BIOMATERIAL WORLD) CLASS DESCRIPTION: BMWorld is a derived class of the parent class World.
 *                                                The BMWorld class manages the model world.
 *                                                It is used to initialize cells, ECM, patches, and chemicals; to destroy
 *                                                agent ArrayChains; to execute each timestep of the model; to sprout agents;
 *                                                to count patches; and to output data.
 */
class BMWorld: public World {
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
    BMWorld(double width = 5, //mm
    		double length = 4, //mm
    		double height = 3, //mm
    		double plength = 0.01 //mm (10 um)
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
    void destroyCell(Cell* &agent);

    /*
     * Description:	Assign a patch type to each patch within bounds of type 
     *
     * Return: void
     * Parameters: void
     */
    void assignPatches(int type, int xmin, int xmax, int ymin, int ymax, int zmin, int zmax);

    /*
     * Description:	Assign a patch type to each patch in the world in row major index manner
     *
     * Return: void
     * Parameters: void
     */
    void initializePatches();

    /*
     * Description:	Helper function update WorldChem TotalChem variables for testing
     *
     * Return: void
     * Parameters: void
     */
    void updateTotalChem(); 

    /*
     * Description:	Initializes chemical concentrations in each patch and initializes the total concentration of each chemical
     *
     * Return: void
     * Parameters: void
     */
    void initializeChemCPU(); 

#ifdef GPU_DIFFUSE
    /*
     * Description:	Initializes chemical concentrations in each patch and initializes the total concentration of each chemical
     *
     * Return: void
     * Parameters: void
     */
    void initializeChemGPU(); 
#endif

    /*
     * Description:	Initializes chemical concentrations in each patch and initializes the total concentration of each chemical
     *
     * Return: void
     * Parameters: void
     */
    void initializeChem(); 

    /*
     * Description:	Initializes all cells to their correct patches
     *
     * Return: void
     * Parameters: void
     */
    void initializeCells();

    /*
     * Description:	Initializes collagen, aggrecan and hyaluronan to their correct patches
     *
     * Return: void
     * Parameters: void
     */
    void initializeECM(); 

    /*
     * Description:	Initializes all damaged patches. Sprouts platelets and fragments ECM on the damaged patches.
     *
     * Return: void
     * Parameters: void
     */
    void initializeDamage();

    /*
     * Description:	Each call to go() simulates 30 minutes or 'real-world' time of the biological model. 
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
     * 							   Pass in false if should be sprouted in tissue (default)
     */
    void sproutAgent(
    		int num,
    		int patchType,
    		int agentType,
    		int xmin, int xmax,
    		int ymin, int ymax,
    		int zmin, int zmax
    		);

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
     *                             Pass in false if should be sprouted in tissue (default)
     */
    void sproutAgentInArea(
    		int num,
    		int patchType,
    		int agentType,
    		int xmin, int xmax,
    		int ymin, int ymax,
    		int zmin, int zmax
    		);

    /*
     * Description:	Function for sprouting cells in the whole world
     *
     * Return: void
     *
     * Parameters: num          -- Number of cells to sprout
     *             patchType    -- Type of patches of sprout on
     *             agentType    -- Type of agent to sprout
     *             bloodOrTiss  -- Pass in true if should be sprouted in blood
     *                             Pass in false if should be sprouted in tissue (default)
     */
    void sproutAgentInWorld(
    		int num,
    		int patchType,
    		int agentType
    		);

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
     * Description:	Converts a length in millimeters to the number of patches
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
     * Description:	Reads user input from config file (default: config.txt) to initialize chemicals, wound, cells.
     *
     * Return: Returns 0 if function proceeded to completion, for testing. Can be removed later.
     *
     * Parameters: void
     */
    int userInput();

//#ifdef MODEL_SCAFFOLD
    /*
     * Description: Calculate Hydrogel initial Elastic modulus, Crosslink density, pore size
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

    /*
     * Description:Degrade Ca-Alg Patches and replace with immature tissue patch (no new ECM produced)
     *
     * Return:
     *
     * Parameters:        numOfPatches  -- number of CaAlg patches to "degrade" and be replaced with patch type tissue                  
     */
    void degradeCaAlg(int numOfPatches);

    /*
     * Description: Print out extra info for debugging purposes
     *
     * Return:
     * Parameters:
     */
    void debugInfo();
//#endif

/****************************************************************
 * OUTPUT SUBROUTINES & VISUALIZATION                           *
 ****************************************************************/

    ///*
    // * Description:	Outputs cell counts and cytokine levels from the current tick to the file "Output/Output_Biomarkers.csv".
    // *              Used for testing.
    // *
    // * Return: void
    // * Parameters: void
    // */
    //void outputWorld_csv();

    /*
     * Description:	Outputs all patch assignments (patch type, agent type, ECM type) to files in output directory.
     *
     * Return: void
     * Parameters: void
     */
    void patchassign_csv();

/****************************************************************
 * STATIC VARIABLES                                             *
 ****************************************************************/

    static unsigned seed;    // Used to generate random numbers
    static bool highTNFdamage;    // Whether there is high TNF damage (which results in ECM fragmentation)
    static float patchpermm;      // The number of patches per millimeter in the world
    static float liveCells;
    static float deadCells;
    static float deletedCells;
    static float prevCells;
    static int initialCaAlg;      // The number of initial tissue patches

    //Crosslinked Ca-Alg Hydrogel Parameters
    // Ca Crosslinker:
    static float Ca_Mw;	 // Molecular Weight 
    // Alginate:
    static float Alg_Mn; // Number Average Molecular Weight (g mol-1)
    static float totalVolumeML; //
    

    /* CALIBRATION Variables */
    static float thresholdTNFdamage;      // The threshold for TNF damage
    static float sproutingFrequency[6];   // The number of hours between agent sprouting sessions
    static float sproutingAmount[14];     // Constants related to the number of agents to sprout
    static float cytokineDecay[6];        // The decay rates of the cytokines
    static float halfLifes_static[6];     // The half lifes of the cytokines in minutes

    /* Calibration Variables */
    static float ElasticMod[7];   // Elastic Modulus of Ca-Alg Hydrogel
    static float XLDensity[2];    // Crosslink density of Ca-Alg Hydrogel
    static float SwellRatio[5];   // Swell Ratio of Ca-Alg Hydrogel
    static float MassLoss[4];     // Mass loss of Ca-Alg Hydrogel
    static float PoreSize[2];     // PoreSize of Ca-Alg Hydrogel


/****************************************************************
 * CONSTANT VARIABLES                                           *
 ****************************************************************/
    double patchlength;    // The length of each patch
  
    // Instance of type to manage chemicals in the world:
    Chemical WorldChem;
    /* Used to allocate the chemicals as floats. It is an array of float*, where each float* points to array of concentrations of a given chemical at different patches. 
     * chemAllocation[chemIndex][patchIndex] can be used as a multidimensional array to access a given chemical concentration using indexing. It is linked to WorldChem. */
    float** chemAllocation;

    float E;    // Elastic Modulus (Pa)
    float pXL;  // Crosslink Density (mmol/mL)
    float Q;    // Swelling Ratio (%)
    float w;    // Mass Loss (%)
    float poreWidth;    // Pore Size (um)

#ifdef GPU_DIFFUSE		// Diffusion using Third buffer
    /* 
    typedef struct CCTX		// convolution context
    {
	int KH;
	int KW;
	int KX;
	int KY;
	int DH;
	int DW;
	int FFTH;
	int FFTW;
    } c_ctx;
    */
    typedef struct CCTX c_ctx;
    float**   h_diffusion_results;
    fComplex** d_kernel_spectrum;
    fComplex** h_dKernel_spectrum;
    c_ctx*     chem_cctx;
#endif

    int typesOfChem;    // The number of different chemicals there are in the world
    vector<float> baselineChem;     // Initial amount of each chemical in the world
    ECM* worldECM;                  // Pointer to array of ECM
    vector<int> initialCells;       // Initial amount of each cell type (agent) in the world

    ArrayChain<Cell*> cells;    // ArrayChain to manage all cell data
    vector<Cell*>* localNewCells[MAX_NUM_THREADS];    // Vector of pointers to local lists of cell pointers to add to global list
    vector<int> initHAcenters;          // Vector of patches which can be centers for sprouting original hyaluronan
    unsigned seeds[MAX_NUM_THREADS];    // Seeds used to generate random numbers for each thread
    float *D;          // Array of diffusion coefficients, gets allocated in userInput()
    int *HalfLifes;    // Arsray of cytokine half-life (seconds), gets allocated in userInput()

    vector<float> tgfLine; // Vector to store TGF values along an x-face line from boundary to center of ABM grid

    int lineY;
    int lineZ;

    float Alg_v, Alg_wv;  // Volume (mL) and final concentration (% w/v) of Alg in Ca-Alg hydrogel
    float Ca_v, Ca_wv;    // Volume (mL) and final concentration (% w/v) of Ca 3400
    float highMW_alg, lowMW_alg; // ratio components of high and low MW kDa in the alginate hydrogel


 protected:
     // --- output file hooks ---
     std::string get_output_filename()                                override;
     void        write_csv_header(std::ofstream& file)               override;
     void        write_data_row(std::ofstream& file,
         std::map<std::string, int>& agent_counts,
         std::map<std::string, float>& env_counts) override;

     // --- auxiliary output hooks ---
     void write_auxiliary_header()                                    override;
     void write_auxiliary_outputs()                                   override;

     // --- agent counting hooks ---
     std::vector<std::string> get_agent_type_names()                   override;
     void count_agent_types(std::map<std::string, int>& agent_counts)   override;

     // --- agent population hooks ---
     int  get_total_agent_count() override;

     // --- environment element counting hooks (e.g. ecm) ---
     std::vector<std::string> get_env_type_names()                    override;
     void count_env(std::map<std::string, float>& env_counts)        override;


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
    //void seedCells(float hours);

    /*
     * Description:	(Stage 1)	Entry point function for diffusing chemicals
     * 							Selects the appropriate diffusion function to apply
     *
     * Return: void
     *
     * Parameters: 
     */
    void diffuseCytokines();

    /*
     * Description:	Helper function for cell execution. Execute all alive chondrocytes.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline runCells();

    /*
     * Description:	(Stage 2)	Execute cell functions for all living cells.
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
     * Description:	Helper function for diffuseCytokines(). 
     *              Update chemical  concentration of all chemicals at each patch over a given time step using Netlogo diffusion procedure.
     *
     * Return: void
     *
     * Parameters: void
     */
    void NetlogoDiffuse();

    /*
     * Description:	Helper function for diffuseCytokines(). 
     *              Update chemical concentration of all chemicals at each patch over a given
     * 				time step using partial differential diffusion equation.
     *
     * Return: void
     *
     * Parameters: coeff  -- Diffusion coefficient (mm^2/min)
     *             dt     -- Time step ( min) 
     *             NOTE: Examples of diffusion coefficients can be found at http://www.math.ubc.ca/~ais/website/status/diffuse.html
     */
    #ifdef MODEL_3D
        // 3D stability condition, dt < dx^2/6*D min = 0.020833
        //assuming dx = dy = dz = 0.015 mm, D = 0.0018 mm^2/min :
        void diffuseChem(int ichem, float dt = 0.02, float coeff = 0.0018);      
    #else
        // 2D stability condition, dt < dx^2/4*D min = 0.03125
        // assuming dx = dy = 0.015 mm, D = 0.0018 mm^2/min
        void diffuseChem(int ichem, float dt = 0.03, float coeff = 0.0018);
    #endif

    
#ifdef GPU_DIFFUSE
    /*
     * Description:	Helper function for diffuseCytokines(). 
     *              Convolution-based chemical diffusion executed on GPU.
     *
     * Return: void
     *
     * Parameters: void
     *             Note:  All parameters are assumed to have been intialized in cctx_t (convolution context) via initializeChemGPU()
     */
    void diffuseChemGPU();
#endif

    /*
     * Description:	(Stage 4a)	Update chemicals to reflect next tick's states
     * 			Differ from updateChem() since this should update in the following manner:
     * 				p<chem> = d<chem> + t<chem>*(1-gamma)
     * 			where gamma is a cytokine specific constant derived from the cytokine's halflife.
     * 				gamma = 1 - 2^(-1/halflife)
     * Return: void
     *
     * Parameters: void
     */
    void updateChem();

    /*
     * Description:	Update chemicals to reflect next tick's states.
     * 			This is called once at the initilzation state.
     * 			Differ from updateChem() since this should update in the following manner:
     * 				p<chem> = d<chem> + p<chem>*(1-gamma)
     * 			where gamma is a cytokine specific constant derived from the cytokine's halflife.
     * 				gamma = 1 - 2^(-1/halflife)
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateChemCPU();
    
    /*
     * Description:	Helper function for ECM updates. Execute updates for ALL ECM managers.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline executeAllECMUpdates();

    /*
     * Description:	Helper function for ECM updates. Execute request resets for ALL ECM managers.
     *
     * Return: void
     *
     * Parameters: void
     */
    void inline executeAllECMResetRequests();

    /*
     * Description:	(Stage 4c)	Update ECM managers to reflect next tick's states
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateECMManagers();

    /*
     * Description:	(Stage 4b)	Update cells to reflect next tick's states
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
    // used only in write_data_row — not exposed as hooks
    float calculate_viability();
    float calculate_pct_differentiated(std::map<std::string, int>& agent_counts);
};
#endif	/* BMWorld_H */