/* 
 * File: Agent.h
 *
 * File Contents: Contains declarations for the Agent class.
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#ifndef AGENT_H
#define	AGENT_H

#include <stdlib.h>
#include <algorithm>
#include <vector>

#include "../common.h"
#include "../enums.h"
//#include "../World/World.h"

class World;
class WHWorld;
class Patch;
class ECM; 

using std::vector;

/*
 * AGENT CLASS DESCRIPTION:       Agent is a parent class of all cell "agents". 
 *                                It is used to move agents on 2D patches, to make agents die, and to make them move to highest chemical concentration.
 */
class Agent {
 public:
    /*
     * Description:	Default agent constructor. 
     *
     * Return: void
     * Parameters: void
     */
    Agent();

    /*
     * Description:	Virtual agent destructor.
     *
     * Return: void
     * Parameters: void
     */
    virtual ~Agent();

    enum agenttype_t {stem, progen, np};            // Enumic type to keep track of the type of agent
    enum chemtype_t {TNF, TGF, IL1beta};            // Enumic type to keep track of the type of chemical. 

    /*
     * Description:	Determines whether an agent is alive or not and updates the agent's properties for the end of the current tick.
     *
     * Return: True if the agent is alive, false otherwise.
     * Parameters: void
     */
    bool isAlive();

    /*
     * Description:	Determines the location of the agent in the x, y, z dimensions of the world.
     *
     * Return: x-coordinate of agent
     * Parameters: void
     */
    int getX();
    int getY();
    int getZ();

    /*
     * Description:	Determines the patch row major index of the agent.
     *
     * Return: patch row major index of agent
     * Parameters: void
     */
    int getIndex();

#ifdef MODEL_SCAFFOLD
    /*
     * Description:	Calculates the speed (patch/tick) cell moves in Ca-Alg hydrogel 
     *
     * Return: void
     * Parameters: void             
     */
    static void calculateMigrationSpeed(int agentType);

    /*
     * Description:	Calculates the rate of duplication per hour in Ca-Alg hydrogel 
     *
     * Return: void
     * Parameters: void            
     */
    static void calculateProliferationRate();

    /*
     * Description:	Calculates the proportion of dead cells 
     *
     * Return: void
     * Parameters: void             
     */
    static void calculateViabilityRate();

    /*
     * Description:	Calculates the amount of collagen, aggrecan and HA produced by cells on Ca-Alg Gel
     *
     * Return: void
     * Parameters: void            
     */
    void calculateECMSynthesisRate(int agentType);
    
    /*
     * Description:	Calls subroutines calculating parameters characterizing cell behavior in Ca-Alg
     *
     * Return: void
     * Parameters: void            
     */
    static void cellCaAlgBehavior();
    
#endif //MODEL_SCAFFOLD

    /*
     * Description:	Rolls a hypothetical dice with 'percent' chance of successful roll
     *
     * Return: True if a successful roll, false otherwise
     *
     * Parameters: percent  -- percentage that dictates the chance of a successful roll
     */
    static bool rollDice(float percent);

    /*
     * Description:	Moves the cell to the target patch if it is valid and available.
     *              Updates cell location and patch availability.
     *
     * Returns: True if the move was successful, false otherwise.
     * 
     * Parameters: dX -- distance to move in the x direction (-: left, +: right)
     *			   dY -- distance to move in the y direction (-: up, +: down)
     *			   dZ -- distance to move in the z direction (-: inwards, +: outwards)
     *	    readIndex -- index into cell's current location array
     *	    NOTE: If the location has been modified in this tick, readIndex should be write_t. Else, readIndex should be read_t
     *
     * Usage: This function should be called after a call to check whether the location is dirty (modified in this tick)
     *        Example call sequence:
     *				    if (isModified(this->index)) move (dX, dY, dZ, write_t);
     *				    else move (dX, dY, dZ, read_t);
     *				    
     *        The check is to make sure that we are reading from the most up-to-date information of the cell. This check is necessary,
     *        since move() can be called multiple times in a tick, thus could potentially modify each attribute more than once. 
     *        The intermediate values need to be kept track of.
     */
    bool move(int dX, int dY, int dZ, int readIndex);

    /*
     * Description:	Determines the mean concentration of chemical of type chemIndex from 27 neighboring patches
     *
     * Return: mean concentration of chemical of type chemIndex from 27 neighboring patches
     *
     * Parameters: chemIndex  -- enumic value of chemical_t for chemical type
     */
    float meanNeighborChem(int chemIndex);

    /*
     * Description:	Determines the number of cells of type cellIndex from 27 neighboring patches.
     *
     * Return: number of cells of type cellIndex from 27 neighboring patches.
     *
     * Parameters: cellIndex  -- enumic value of agent_t for cell type
     */
    int countNeighborCells(int cellIndex);

    /*
     * Description:	Determines the number of ECM proteins of type ECMIndex from 27 neighboring patches.
     *
     * Return: number of ECM proteins of type ECMIndex from 27 neighboring patches.
     *
     * Parameters: ECMIndex  -- enumic value of agent_t for ECM protein type
     */
    int countNeighborECM(int ECMIndex);

    /*
     * Description:	Moves to neighboring patch with highest concentration of chemical of type chemIndex.
     *
     * Return: True if moved successfully, false otherwise.
     *
     * Parameters: chemIndex  -- enumic value of chemical_t for chemical type
     */
    bool moveToHighestChem(int chemIndex);

    /*
     * Description:	Move to a random neighboring patch while respecting rules dictating movement between tissue types.
     *
     * Return: True if moved successfully, false otherwise.
     *
     * Parameters: chemIndex  -- enumic value of chemical_t for chemical type
     */
    void wiggle();

    /*
     * Description:	Updates class members for the end of the tick.
     *
     * Return: void
     * Parameters: void
     */
    void updateAgent();

    /*
     * Description:	Virtual function for cell function
     *
     * Return: void
     * Parameters: void
     */
    virtual void cellFunction();

    /*
     * Description:	Virtual function for agent death
     *
     * Return: void
     * Parameters: void
     */
    virtual void die() = 0;

    /* 
     * Description: Virtual function for copying and initializing a new agent
     *
     * Return: void
     *
     * Parameters: original  -- Agent to be copied
     *             dx        -- Difference in x-coordinate between current and new agents.
     *             dy        -- Difference in y-coordinate between current and new agents.
     *             dz        -- Difference in z-coordinate between current and new agents.
     *                          NOTE: dz = 0 because it is only 2D for now.
     */
    virtual void copyAndInitialize(Agent* original, int dx, int dy, int dz = 0);

    /*
    * Description:	Hatches a new cell on 'number' unoccupied neighbors.
    *              Does not update numOfCells; this must be done elsewhere.
    *
    * Return: void
    *
    * Parameters: number  -- Number of new cells to hatch
    */
    virtual void hatchnewcell(int number, int agentType);

/* -------------------------------------------------------------------------- */
/*                              PUBLIC VARIABLES                              */
/* -------------------------------------------------------------------------- */
    int life[2];        // Number of lives remaining at the beginning and end of each tick
    bool activate[2];   // Whether agent is activated or not at the beginning and end of each tick
    int color[2];       // Agent's color at the beginning and end of each tick
    float size[2];      // Agent's size at the beginning and end of each tick
    int type[2];        // Agent's type at the beginning and end of each tick

/* -------------------------------------------------------------------------- */
/*                              STATIC VARIABLES                              */
/* -------------------------------------------------------------------------- */
    static WHWorld* agentWorldPtr;  // Pointer from an agent to a WHWorld  
    static Patch* agentPatchPtr;    // Pointer from an agent to a Patch
    static ECM* agentECMPtr;        // Pointer from an agent to an ECM
    static int nx, ny, nz;          // Number of patches in x,y,z dimensions of the world
    static int dX[27], dY[27], dZ[27];   // Difference in x-,y-,z-coordinates to neighboring patches

    #ifdef MODEL_3D
        static int neighbor[27]; // Array of neighboring patches in 3D model
    #else 
        static int neighbor[8];  // Array of neighboring patches in 2D model
    #endif

/* -------- Parameters related to chondrocyte behavior in Ca-Alg Gel -------- */
#ifdef MODEL_SCAFFOLD
    static bool CaAlgFlag;          // Flag indicating if static parameters have been set             
    static float proliferationRate; // Change in population (% of initial population) over 1 hour
    static float viabilityRate;     // Viability Rate (%) of cells 
    static float HASynthRate;       // Amount of HA synthesized in Ca-Alg(10^-4 ug)
#endif

/* -------------------------- Calibration Variables ------------------------- */
#ifdef MODEL_SCAFFOLD
    static float CaAlgProlif[5];     // Parameters invloved in Chondrocyte proliferation in CaAlg Gel
    static float CaAlgViability[3];  // Parameters invloved in Chondrocyte viability in CaAlg Gel
    static float HASynth[3];         // Parameters invloved in HA synthesis in CaAlg Gel
#endif

/* --------------------------- PROTECTED VARIABLES -------------------------- */
 protected:
    int ix[2],iy[2],iz[2];  // Agent position in x,y,z dimensions at the beginning and end of each tick
    int index[2];           // Patch row major index for agent at the beginning and end of each tick
    bool alive[2];          // Life status of agent at the beginning and end of each tick
};

#endif	/* AGENT_H */