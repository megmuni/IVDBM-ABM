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
#include "../../World/Usr_World/woundHealingWorld.h"

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
     * Description:	Moves a cell along its preferred chemical gradient.
     *
     * Return: void
     * Parameters: void
     */
    void cellSniff();

    /*
     * Description:	Performs cell death. Updates the cell class members. 
     *              Does not update numOfCells; this must be done elsewhere.
     *
     * Return: void
     * Parameters: void
     */
    void die();						

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
     * Description:	Sprouts new hyaluronan on one of the activated chondrocyte's damaged neighbor patches.
     *
     * Return: void
     *
     * Parameters: meanTNF  -- Average TNF concentration of the activated chondrocyte's neighbors
     * 				     meanTGF  -- Average TGF concentration of the activated chondrocyte's neighbors
     * 				     meanIL1  -- Average IL1 concentration of the activated chondrocyte's neighbors
     */
    //void makeHyaluronan(float meanTNF, float meanTGF, float meanIL1); /*NOTE: not used in stem cell biomaterial ABM

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
    static float cytokineSynthesis[10];   // Parameters involved in synthesis of TNF, TGF, IL1beta by cells
    static float activation[5];           // Parameters involved in cell activation and deactivation
    static float ECMsynthesis[12];        // Parameters involved in ECM synthesis
    static float proliferation[6];        // Parameters invloved in cell proliferation

  protected:
      /*
      * Description:	Determines whether an agent is still proliferative (i.e. has it reached the max number of doublings).
      * Added by MM, 2025.
      * 
      * Return: True if the agent is proliferative, false otherwise.
      * Parameters: agentType
      */
      virtual bool isProliferative();
      virtual int get_max_doublings(); // gets maximum number of cell divisions depending on cell type
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
  
    /*
       * Description:	Performs biological function of a stem cell.
       *
       * Return: void
       * Parameters: void
       */                                                                                                
    void stem_cellFunction(); 
  
    /*
       * Description:	Differentiates the stem cell to the next stage (progenitor).
       *              Does not update numOfStem; this must be done elsewhere.
       *
       * Return: void
       *
       * Parameters: void
       */
    void differentiateStem(int number, int agentType);
  
  /* -------------------------------------------------------------------------- */
  /*                              STATIC VARIABLES                              */
  /* -------------------------------------------------------------------------- */
  static int numOfStem; // Keeps track of the quantity of living stem cells
  static float migrationSpeed;    // Speed (patch/tick) stem cells move in world

  static float collagenSynthRate; // Amount of collagen synthesized in Ca-Alg(10^-4 ug)
  static float aggrecanSynthRate; // Amount of aggrecan synthesized in Ca-Alg(10^-4 ug)
  static float apoptosisChance;
  static float divisionNum; // number of cell divisions the cell has undertaken

  /* -------------------------- Calibration variables ------------------------- */
  static float CaAlgMigration[2];  // Parameters invloved in stem cell migration speed in CaAlg Gel
  static float cytokineSynthesis[3]; // Parameters involved in cytokine synthesis by stem cells (baseline rates)
  static float CollagenSynth[1];   // Parameters invloved in collagen synthesis in CaAlg Gel
  static float AggrecanSynth[1];   // Parameters invloved in aggrecan synthesis in CaAlg Gel
  static float ECMsynthesis[4]; // Parameters involved in ECM synthesis (baseline rates, hours between synth)
  static float proliferation[5]; // Parameters involved in stem cell proliferation (coefficients for probabilistic differentiation
  static float differentiation[5]; // Parameters involved in stem cell differentiation
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

    /*
       * Description:	Performs biological function of an NP progenitor cell.
       *
       * Return: void
       * Parameters: void
       */                                                                                                
    void progen_cellFunction(); 
  
    /*
       * Description:	Differentiates the progenitor cell to the next stage (progenitor).
       *              Does not update numOfProgen; this must be done elsewhere.
       *
       * Return: void
       *
       * Parameters: void
       */
    void differentiateProgen(int number, int agentType);
  
  /* -------------------------------------------------------------------------- */
  /*                              STATIC VARIABLES                              */
  /* -------------------------------------------------------------------------- */
  static int numOfProgen; // Keeps track of the quantity of living progenitor cells
  static float migrationSpeed;    // Speed (patch/tick) pre-NP cells move in world
  static float apoptosisChance;
  static float divisionNum; // number of cell divisions the cell has undertaken

  //static float collagenSynthRate; // Amount of collagen synthesized in Ca-Alg(10^-4 ug) // may not be applicable to Pre-NP
  static float aggrecanSynthRate; // Amount of aggrecan synthesized in Ca-Alg(10^-4 ug)

  /* -------------------------- Calibration variables ------------------------- */
  static float CaAlgMigration[2];  // Parameters invloved in pre-NP cell migration speed in CaAlg Gel
  static float cytokineSynthesis[3]; // Parameters involved in cytokine synthesis by pre-NP cells (baseline rates)
  static float AggrecanSynth[1]; // Parameters involved in ECM synthesis (baseline rates, hours between synth)
  static float proliferation[1]; // Parameters involved in pre-NP cell proliferation (coefficients for probabilistic differentiation
  static float differentiation[3]; // Parameters involved in pre-NP cell differentiation

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

    /*
       * Description:	Performs biological function of an NP cell.
       *
       * Return: void
       * Parameters: void
       */                                                                                                
    void NP_cellFunction();
  
  /* -------------------------------------------------------------------------- */
  /*                              STATIC VARIABLES                              */
  /* -------------------------------------------------------------------------- */
  static int numOfNP; // Keeps track of the quantity of living NP cells
  static float migrationSpeed;    // Speed (patch/tick) NP cells move in world
  static float divisionNum; // number of cell divisions the cell has undertaken

  static float collagenSynthRate; // Amount of collagen synthesized in Ca-Alg(10^-4 ug)
  static float aggrecanSynthRate; // Amount of aggrecan synthesized in Ca-Alg(10^-4 ug)

  /* -------------------------- Calibration variables ------------------------- */
  static float CaAlgMigration[2];  // Parameters invloved in NP cell migration speed in CaAlg Gel
  static float CollagenSynth[3];   // Parameters invloved in collagen synthesis in CaAlg Gel
  static float AggrecanSynth[3];   // Parameters invloved in aggrecan synthesis in CaAlg Gel
};

#endif