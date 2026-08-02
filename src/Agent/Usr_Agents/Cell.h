/* 
 * File: Cell.h
 *
 * File Contents: Contains declarations for the Cell class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *               Meghana Munipalle
 */

#ifndef CELL_H
#define	CELL_H

#include "../Agent.h"
#include "../../Patch/Patch.h"
#include "../../World/Usr_World/biomaterialWorld.h"
#include "../biology_parameters_config.h"

class ECM; 

#include <stdlib.h>
#include <vector>
#include <cmath>
#include <omp.h>
//#include <memory>

using namespace std;

/*
 * CELL CLASS DESCRIPTION:    Cell is a derived class of the parent class Agent. 
 *                                   It manages all cell agents, including stem, progenitor, and NP.
 *                                   It is used to initialize a cell, carry out its biological function,
 *                                   activate it, deactivate it, make it move, and make it die.
 */
class Cell: public Agent {
  public:
    /*
     * Description:	Default cell constructor. 
     *
     * Return: void
     * Parameters: void
     */
    Cell();

    /*
     * Description:	Cell constructor. Initializes cell attributes.
     *
     * Return: void
     * Parameters: patchPtr  -- Pointer to patch on which cell will reside.
     *                          NOTE: The pointer cannot be NULL.
     */ 
    Cell(Patch* patchPtr);

    /*
     * Description:	Cell constructor. Initializes cell class members.
     *
     * Return: void
     * 
     * Parameters: x  -- Position of cell in x dimension
     *             y  -- Position of cell in y dimension
     *             z  -- Position of cell in z dimension
     */
    Cell(int x, int y, int z); 

    /*
     * Description:	Cell destructor.
     *
     * Return: void
     * Parameters: void
     */
    ~Cell(); 

    /*
     * Description:	Performs biological function of a cell.
     *
     * Return: void
     * Parameters: void
     */
    void cellFunction();

    /*
     * Description:	Performs biological function of an unactivated chondrocyte.
     *
     * Return: void
     * Parameters: void
     */                                                                                                
    /*void chond_cellFunction();                                                                                                                        

    /*
     * Description:	Performs biological function of an activated chondrocyte.
     *
     * Return: void
     * Parameters: void
     */
    /*void achond_cellFunction();                         

    /*
     * Description:	Moves a cell along its preferred chemical gradient. Template method.
     *
     * Return: void
     * Parameters: void
     */
    virtual void cellSniff() final;

    /*
     * Description:	Performs cell death. Updates the cell class members. 
     *              Does not update numOfCells; this must be done elsewhere.
     *
     * Return: void
     * Parameters: void
     */
    void die();					

    /*
     * Description:	Spefically carries out cell apoptosis based on apoptosis chances for cell types
     * 
     *
     * Return: void
     * Parameters: void
     */
    virtual void apoptose() final;

    /*
     * Description:	Activates an unactivated chondrocyte. Updates the chondrocyte class members.
     *
     * Return: void
     * Parameters: void
     */
    /*void chondActivation();

    /*
     * Description:	Deactivates an activated MSC. Updates the cell class members.
     *
     * Return: void
     * Parameters: void
     */
    //void cellDeactivation();

    /* 
     * Description: Copies the location of 'original' agent and initializes a new chondrocyte at a distance away determined by dx, dy, dz.
     *              NOTE: Target patch at distance of dx,dy,dz must be unoccupied for proper functionality.
     *
     * Return: void
     *
     * Parameters: original  -- Agent to be copied
     *             dx        -- Difference in x-coordinate of the cell's location relative to original's location.
     *             dy        -- Difference in y-coordinate of the cell's location relative to original's location.
     *             dz        -- Difference in z-coordinate of the cell's location relative to original's location.
     *                          NOTE: dz = 0 because it is only 2D for now.
     */
    void copyAndInitialize(Agent* original, int dx, int dy, int dz = 0);

    /*
     * Description: template method for cell proliferation
     * 
     */
    virtual void proliferate() final;

    /*
     * Description: template method for cell differentiation. Called once per simulation
     * tick for each cell. Checks all conditions before attempting cell division, then creates
     * a daughter cell and increments that cell's number of divisions.
     *
     */
    virtual void differentiate() final;

    /*
     * Description: template method for ECM synthesis by cells
     *
     */
    virtual void ecm_synthesis() final;

    /*
     * Description:	Sprouts original collagen on one of the activated cell's damaged neighbor patches.
     *
     * Return: void
     *
     * Parameters: meanTGF   -- Average TGF concentration of the cell's neighbors
     * 				     meanIL1   -- Average IL1 concentration of the cell's neighbors
     * 				     countnHA	 -- Number of new hyaluronan on the cell's neighbors
     * 				     countfHA	 -- Number of fragment hyaluronan on the cell's neighbors
     */
    void makeOCollagen(float meanTGF, float meanIL1); //, int countnHA, int countfHA

    /*
     * Description:	Sprouts original aggrecan on one of the cell's damaged neighbor patches.
     *
     * Return: void
     *
     * Parameters: meanTNF  -- Average TNF concentration of the cell's neighbors
     * 				     meanTGF  -- Average TGF concentration of the cell's neighbors
     * 				     meanIL1  -- Average IL1 concentration of the cell's neighbors
     */
    void makeOAggrecan(float meanTNF, float meanTGF, float meanIL1);

    /*
     * Description: template method for cytokine synthesis by cells
     *
     */
    virtual void cytokine_synthesis() final;

    /*
    * Description:	Hatches a new cell on 'number' unoccupied neighbors.
    *              Does not update numOfCells; this must be done elsewhere.
    *
    * Return: void
    *
    * Parameters: number    -- Number of new cells to hatch
    *             agentType -- Type of agent/cell to hatch
    *             here      -- Whether to hatch on neighboring patches (default 0) or current patch (1)
    */
    virtual void hatchnewcell(int number, int agentType, int here = 0);


    /* -------------------------------------------------------------------------- */
    /*                              STATIC VARIABLES                              */
    /* -------------------------------------------------------------------------- */
    static int numOfCells;  // Keeps track of the quantitiy of living cells.

    /* -------------------------- Calibration variables ------------------------- */
    enum ProliferationIdx {
      PROLIFERATION_HOURS_BETWEEN = 0,
      PROLIFERATION_TGF_THRESHOLD,
      PROLIFERATION_LOG_SCALE,
      PROLIFERATION_LOG_OFFSET,
      PROLIFERATION_COUNT
    };
    static_assert(sizeof(CellProliferationParams) / sizeof(double) == PROLIFERATION_COUNT,
                  "CellProliferationParams field count must match Cell::ProliferationIdx");
    enum CytokineSynthesisIdx {
      CYTOKINE_TGF_BASELINE = 0,
      CYTOKINE_TGF_FEEDBACK_TGF,
      CYTOKINE_TGF_FEEDBACK_IL1BETA,
      CYTOKINE_TGF_FEEDBACK_TNF,
      CYTOKINE_TNF_BASELINE,
      CYTOKINE_TNF_FEEDBACK_IL1BETA,
      CYTOKINE_TNF_FEEDBACK_TGF_DENOM,
      CYTOKINE_IL1BETA_BASELINE,
      CYTOKINE_IL1BETA_FEEDBACK_TNF,
      CYTOKINE_IL1BETA_FEEDBACK_TGF_DENOM,
      CYTOKINE_SYNTHESIS_COUNT
    };
    static_assert(sizeof(CellCytokineSynthesisParams) / sizeof(double) == CYTOKINE_SYNTHESIS_COUNT,
                  "CellCytokineSynthesisParams field count must match Cell::CytokineSynthesisIdx");
    static float cytokineSynthesis[CYTOKINE_SYNTHESIS_COUNT];   // Parameters involved in synthesis of TNF, TGF, IL1beta by cells
    static float proliferation[PROLIFERATION_COUNT];        // Parameters invloved in cell proliferation

  protected:

      // Proliferation-related hook functions
      /*
      * Description:	Determines whether an agent is still proliferative (i.e. has it reached the max number of doublings).
      * Added by MM, 2025.
      * 
      * Return: True if the agent is proliferative, false otherwise.
      * Parameters: agentType
      */
      virtual bool isProliferative(); // checks if a cell can proliferate based on its current number of doublings
      virtual float get_prolif_prob(float meanTGF,
          float meanIL1,
          float meanTNF); // calculates probability of dividing based on local chemical
      virtual int get_max_doublings(); // gets maximum number of cell divisions depending on cell type

      // Differentiation-related hook functions
      virtual int get_daughter_type(); // determines what cell type is produced
      virtual float get_diff_prob(float meanTGF,
          float meanIL1,
          float meanTNF); // calculates probability of differentiating based on local chemical

      // ECM synthesis-related hook functions
      virtual void calculate_ecm_synth_rates(float meanTGF, float meanIL1, float meanTNF, float patchesVolume); // calculates the ECM synthesis rates for each cell
      virtual void create_ecm(float meanTGF, float meanIL1, float meanTNF); // actually creates/produces the ECM on the patch

      // Cytokine-related hook functions
      virtual void create_cytokines(float patchTGF, float patchIL1beta, float patchTNF); // actually creates/produces the cytokines on the patch

      // Movement-related hook functions
      virtual float get_migration_speed(); // calculates the migration speed for each cell type and passes it to the [celltype]::migrationSpeed variable
      virtual bool can_tgf_excite() { return false; } // used only for NP cell to check if TGF can excite the cell into moving

      // Death-related hook functions
      virtual float get_apoptosis_chance(); // calculates the chance of cell death 

      // OCR-related hook functions
      virtual float get_OCR(); // gets the cell-specific OCR 
};

/*
 * STEM CLASS DESCRIPTION:    Stem is a derived class of the parent class Cell. 
 *                                   It manages stem cells and is used to make stem cells differentiate.                         
 */
class Stem: public Cell {
  public:

    /*
     * Description:	Default Stem constructor.
     *
     * Return: void
     * Parameters: void
     */
    Stem();

    /*
     * Description:	Stem constructor. Initializes stem cell attributes.
     *
     * Return: void
     * Parameters: patchPtr  -- Pointer to patch on which stem cell will reside.
     *                          NOTE: The pointer cannot be NULL.
     */
    Stem(Patch* patchPtr);


    /*
     * Description:	Stem constructor. Initializes stem class members.
     *
     * Return: void
     *
     * Parameters: x  -- Position of cell in x dimension
     *             y  -- Position of cell in y dimension
     *             z  -- Position of cell in z dimension
     */
    Stem(int x, int y, int z);

    /*
     * Description:	Stem destructor.
     *
     * Return: void
     * Parameters: void
     */
    ~Stem();
  
  /* -------------------------------------------------------------------------- */
  /*                              STATIC VARIABLES                              */
  /* -------------------------------------------------------------------------- */
  static int numOfStem; // Keeps track of the quantity of living stem cells
  static float migrationSpeed;    // Speed (patch/tick) stem cells move in world

  static float OCR; // Stem cell oxygen consumption rate in fmol/h/cell
  static float collagenSynthRate; // Amount of collagen synthesized in Ca-Alg(10^-4 ug)
  static float aggrecanSynthRate; // Amount of aggrecan synthesized in Ca-Alg(10^-4 ug)
  static float apoptosisChance;
  static float divisionNum; // number of cell divisions the cell has undertaken

  /* -------------------------- Calibration variables ------------------------- */
  enum MigrationIdx { MIGRATION_ELASTICITY_EFFECT = 0, MIGRATION_BASELINE_SPEED, MIGRATION_COUNT };
  static_assert(sizeof(MigrationParams) / sizeof(double) == MIGRATION_COUNT,
                "MigrationParams field count must match Stem::MigrationIdx");
  static float CaAlgMigration[MIGRATION_COUNT];  // Parameters invloved in stem cell migration speed in CaAlg Gel

  enum CytokineSynthesisIdx { CYTOKINE_TGF_BASELINE = 0, CYTOKINE_TNF_BASELINE, CYTOKINE_IL1BETA_BASELINE, CYTOKINE_SYNTHESIS_COUNT };
  static_assert(sizeof(StemCytokineSynthesisParams) / sizeof(double) == CYTOKINE_SYNTHESIS_COUNT,
                "StemCytokineSynthesisParams field count must match Stem::CytokineSynthesisIdx");
  static float cytokineSynthesis[CYTOKINE_SYNTHESIS_COUNT]; // Parameters involved in cytokine synthesis by stem cells (baseline rates)

  static_assert(sizeof(StemCollagenSynthesisParams) / sizeof(double) == 1,
                "StemCollagenSynthesisParams field count must match Stem::CollagenSynth size");
  static float CollagenSynth[1];   // Parameters invloved in collagen synthesis in CaAlg Gel
  static_assert(sizeof(StemAggrecanSynthesisParams) / sizeof(double) == 1,
                "StemAggrecanSynthesisParams field count must match Stem::AggrecanSynth size");
  static float AggrecanSynth[1];   // Parameters invloved in aggrecan synthesis in CaAlg Gel

  enum ProliferationIdx {
    PROLIFERATION_TGF_THRESHOLD = 0,
    PROLIFERATION_TNF_EFFECT,
    PROLIFERATION_IL1BETA_EFFECT,
    PROLIFERATION_ELASTICITY_EFFECT,
    PROLIFERATION_COUNT
  };
  static_assert(sizeof(StemProliferationParams) / sizeof(double) == PROLIFERATION_COUNT,
                "StemProliferationParams field count must match Stem::ProliferationIdx");
  static float proliferation[PROLIFERATION_COUNT]; // Parameters involved in stem cell proliferation (coefficients for probabilistic differentiation

  enum DifferentiationIdx {
    DIFFERENTIATION_ASYMMETRIC_PROBABILITY = 0,
    DIFFERENTIATION_BASELINE_PROBABILITY,
    DIFFERENTIATION_TGF_EFFECT,
    DIFFERENTIATION_HOURS_BETWEEN_ATTEMPTS,
    DIFFERENTIATION_COUNT
  };
  static_assert(sizeof(StemDifferentiationParams) / sizeof(double) == DIFFERENTIATION_COUNT,
                "StemDifferentiationParams field count must match Stem::DifferentiationIdx");
  static float differentiation[DIFFERENTIATION_COUNT]; // Parameters involved in stem cell differentiation

protected:
    int get_max_doublings() override;
    float get_prolif_prob(float meanTGF,
        float meanIL1,
        float meanTNF) override;

    // Differentiation-related hook functions
    int get_daughter_type() override;
    float get_diff_prob(float meanTGF,
        float meanIL1,
        float meanTNF) override;

    // ECM synthesis-related hook functions
    void calculate_ecm_synth_rates(float meanTGF, float meanIL1, float meanTNF, float patchesVolume) override;
    void create_ecm(float meanTGF, float meanIL1, float meanTNF) override;

    // Cytokine-related hook functions
    void create_cytokines(float patchTGF, float patchIL1beta, float patchTNF) override;

    // Movement-related hook functions
    float get_migration_speed() override;
    
    // Death-related hook functions
    float get_apoptosis_chance() override;

    // OCR-related hook functions
    float get_OCR() override;

};

/*
 * PROGEN CLASS DESCRIPTION:    Progen is a derived class of the parent class Cell. 
 *                                   It manages NP progenitor cells and is used to make them differentiate.                         
 */
class Progen: public Cell {
  public:
    /*
     * Description:	Default Pre-NP constructor.
     *
     * Return: void
     * Parameters: void
     */
    Progen();

    /*
       * Description:	Pre-NP constructor. Initializes pre-NP cell attributes.
       *
       * Return: void
       * Parameters: patchPtr  -- Pointer to patch on which pre-NP cell will reside.
       *                          NOTE: The pointer cannot be NULL.
       */
    Progen(Patch* patchPtr);

    /*
     * Description:	Pre-NP constructor. Initializes progen class members.
     *
     * Return: void
     *
     * Parameters: x  -- Position of cell in x dimension
     *             y  -- Position of cell in y dimension
     *             z  -- Position of cell in z dimension
     */
    Progen(int x, int y, int z);
    
    /*
     * Description:	Pre-NP destructor.
     *
     * Return: void
     * Parameters: void
     */
    ~Progen();
  
  /* -------------------------------------------------------------------------- */
  /*                              STATIC VARIABLES                              */
  /* -------------------------------------------------------------------------- */
  static int numOfProgen; // Keeps track of the quantity of living progenitor cells
  static float migrationSpeed;    // Speed (patch/tick) pre-NP cells move in world
  static float apoptosisChance;
  static float divisionNum; // number of cell divisions the cell has undertaken
  static float OCR; // Pre-NP cell oxygen consumption rate in fmol/h/cell
  //static float collagenSynthRate; // Amount of collagen synthesized in Ca-Alg(10^-4 ug) // may not be applicable to Pre-NP
  static float aggrecanSynthRate; // Amount of aggrecan synthesized in Ca-Alg(10^-4 ug)

  /* -------------------------- Calibration variables ------------------------- */
  enum MigrationIdx { MIGRATION_ELASTICITY_EFFECT = 0, MIGRATION_BASELINE_SPEED, MIGRATION_COUNT };
  static_assert(sizeof(MigrationParams) / sizeof(double) == MIGRATION_COUNT,
                "MigrationParams field count must match Progen::MigrationIdx");
  static float CaAlgMigration[MIGRATION_COUNT];  // Parameters invloved in pre-NP cell migration speed in CaAlg Gel

  enum CytokineSynthesisIdx { CYTOKINE_TGF_BASELINE = 0, CYTOKINE_TNF_BASELINE, CYTOKINE_IL1BETA_BASELINE, CYTOKINE_SYNTHESIS_COUNT };
  static_assert(sizeof(ProgenCytokineSynthesisParams) / sizeof(double) == CYTOKINE_SYNTHESIS_COUNT,
                "ProgenCytokineSynthesisParams field count must match Progen::CytokineSynthesisIdx");
  static float cytokineSynthesis[CYTOKINE_SYNTHESIS_COUNT]; // Parameters involved in cytokine synthesis by pre-NP cells (baseline rates)

  static_assert(sizeof(ProgenAggrecanSynthesisParams) / sizeof(double) == 1,
                "ProgenAggrecanSynthesisParams field count must match Progen::AggrecanSynth size");
  static float AggrecanSynth[1]; // Parameters involved in ECM synthesis (baseline rates, hours between synth)
  //static float proliferation[1]; // Parameters involved in pre-NP cell proliferation. Values are the same as Cell/Stem; can just use the equivalent params defined in Stem
  //static float differentiation[3]; // Parameters involved in pre-NP cell differentiation. Values are the same as Stem; can just use the equivalent params defined in Stem

protected:
    int get_max_doublings() override;
    float get_prolif_prob(float meanTGF,
        float meanIL1,
        float meanTNF) override;

    // Differentiation-related hook functions
    int get_daughter_type() override;
    float get_diff_prob(float meanTGF,
        float meanIL1,
        float meanTNF) override;

    // ECM synthesis-related hook functions
    void calculate_ecm_synth_rates(float meanTGF, float meanIL1, float meanTNF, float patchesVolume) override;
    void create_ecm(float meanTGF, float meanIL1, float meanTNF) override;

    // Cytokine-related hook functions
    void create_cytokines(float patchTGF, float patchIL1beta, float patchTNF) override;

    // Movement-related hook functions
    float get_migration_speed() override;

    // Death-related hook functions
    float get_apoptosis_chance() override;

    // OCR-related hook functions
    float get_OCR() override;
};

/*
 * NP CLASS DESCRIPTION:    NP is a derived class of the parent class Cell. 
 *                                   It defines NP cells.                         
 */
class NP: public Cell {
  public:
    /*
     * Description:	Default NP constructor.
     *
     * Return: void
     * Parameters: void
     */
    NP();

    /*
       * Description:	NP constructor. Initializes NP cell attributes.
       *
       * Return: void
       * Parameters: patchPtr  -- Pointer to patch on which NP cell will reside.
       *                          NOTE: The pointer cannot be NULL.
       */
    NP(Patch* patchPtr);

    /*
       * Description:	NP constructor. Initializes NP class members.
       *
       * Return: void
       *
       * Parameters: x  -- Position of cell in x dimension
       *             y  -- Position of cell in y dimension
       *             z  -- Position of cell in z dimension
       */
    NP(int x, int y, int z);

    /*
     * Description:	NP destructor.
     *
     * Return: void
     * Parameters: void
     */
    ~NP();
  
  /* -------------------------------------------------------------------------- */
  /*                              STATIC VARIABLES                              */
  /* -------------------------------------------------------------------------- */
  static int numOfNP; // Keeps track of the quantity of living NP cells
  static float migrationSpeed;    // Speed (patch/tick) NP cells move in world
  static float apoptosisChance;
  static float divisionNum; // number of cell divisions the cell has undertaken

  static float OCR; // NP cell oxygen consumption rate in fmol/h/cell
  static float collagenSynthRate; // Amount of collagen synthesized in Ca-Alg(10^-4 ug)
  static float aggrecanSynthRate; // Amount of aggrecan synthesized in Ca-Alg(10^-4 ug)

  /* -------------------------- Calibration variables ------------------------- */
  enum MigrationIdx { MIGRATION_ELASTICITY_EFFECT = 0, MIGRATION_BASELINE_SPEED, MIGRATION_COUNT };
  static_assert(sizeof(MigrationParams) / sizeof(double) == MIGRATION_COUNT,
                "MigrationParams field count must match NP::MigrationIdx");
  static float CaAlgMigration[MIGRATION_COUNT];  // Parameters invloved in NP cell migration speed in CaAlg Gel

  enum CollagenSynthIdx { COLLAGEN_SCALING_FACTOR = 0, COLLAGEN_TIME_EFFECT, COLLAGEN_BASELINE_RATE, COLLAGEN_SYNTH_COUNT };
  static_assert(sizeof(NpCollagenSynthesisParams) / sizeof(double) == COLLAGEN_SYNTH_COUNT,
                "NpCollagenSynthesisParams field count must match NP::CollagenSynthIdx");
  static float CollagenSynth[COLLAGEN_SYNTH_COUNT];   // Parameters invloved in collagen synthesis in CaAlg Gel

  enum AggrecanSynthIdx { AGGRECAN_SCALING_FACTOR = 0, AGGRECAN_TIME_EFFECT, AGGRECAN_BASELINE_RATE, AGGRECAN_SYNTH_COUNT };
  static_assert(sizeof(NpAggrecanSynthesisParams) / sizeof(double) == AGGRECAN_SYNTH_COUNT,
                "NpAggrecanSynthesisParams field count must match NP::AggrecanSynthIdx");
  static float AggrecanSynth[AGGRECAN_SYNTH_COUNT];   // Parameters invloved in aggrecan synthesis in CaAlg Gel

protected:
    int get_max_doublings() override;
    float get_prolif_prob(float meanTGF,
        float meanIL1,
        float meanTNF) override;

    // ECM synthesis-related hook functions
    void calculate_ecm_synth_rates(float meanTGF, float meanIL1, float meanTNF, float patchesVolume) override;
    void create_ecm(float meanTGF, float meanIL1, float meanTNF) override;

    // Cytokine-related hook functions
    void create_cytokines(float patchTGF, float patchIL1beta, float patchTNF) override;

    // Movement-related hook functions
    float get_migration_speed() override;
    bool can_tgf_excite() override;

    // Death-related hook functions
    float get_apoptosis_chance() override;

    // OCR-related hook functions
    float get_OCR() override;
};

#endif