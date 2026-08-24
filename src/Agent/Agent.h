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
#define AGENT_H

#include <algorithm>
#include <stdlib.h>
#include <vector>

#include "../Chemistry/species_id.h"
#include "../World/World.h"
#include "../common.h"
#include "../enums.h"
#include "../Utilities/rng.h"
// #include "../ArrayChain/ArrayChain.h"

class World;
class BMWorld;
class ChemicalEnvironment;
class Patch;
class ECM;

using std::vector;

/*
 * AGENT CLASS DESCRIPTION:       Agent is a parent class of all cell "agents".
 *                                It is used to move agents on 2D patches, to
 * make agents die, and to make them move to highest chemical concentration.
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

  enum agenttype_t {
    stem,
    progen,
    np
  }; // Enumic type to keep track of the type of agent

  /*
   * Description:	Determines whether an agent is alive or not and updates
   * the agent's properties for the end of the current tick.
   *
   * Return: True if the agent is alive, false otherwise.
   * Parameters: void
   */
  bool isAlive();

  /*
   * Description:	Determines whether an agent is actually naturally dead
   * (this->realDeath) for counting purposes. Added by MM, 2025.
   *
   * Return: True if the agent is alive, false otherwise.
   * Parameters: void
   */
  bool isRealDead();

  /*
   * Description:	Determines the location of the agent in the x, y, z
   * dimensions of the world.
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

  /*
   * Description:	Rolls a hypothetical dice with 'percent' chance of
   * successful roll
   *
   * Return: True if a successful roll, false otherwise
   *
   * Parameters: percent  -- percentage that dictates the chance of a successful
   * roll
   */
  bool rollDice(float percent);

  /*
   * Description:	Moves the cell to the target patch if it is valid and
   * available. Updates cell location and patch availability.
   *
   * Returns: True if the move was successful, false otherwise.
   *
   * Parameters: dX -- distance to move in the x direction (-: left, +: right)
   *			   dY -- distance to move in the y direction (-: up, +:
   * down) dZ -- distance to move in the z direction (-: inwards, +: outwards)
   *	    readIndex -- index into cell's current location array
   *	    NOTE: If the location has been modified in this tick, readIndex
   * should be write_t. Else, readIndex should be read_t
   *
   * Usage: This function should be called after a call to check whether the
   * location is dirty (modified in this tick) Example call sequence: if
   * (isModified(this->index)) move (dX, dY, dZ, write_t); else move (dX, dY,
   * dZ, read_t);
   *
   *        The check is to make sure that we are reading from the most
   * up-to-date information of the cell. This check is necessary, since move()
   * can be called multiple times in a tick, thus could potentially modify each
   * attribute more than once. The intermediate values need to be kept track of.
   */
  bool move(int dX, int dY, int dZ, int readIndex);

  /**
   * Mean patch concentration of a diffusing species over the 3x3x3 neighbor
   * stencil (SpeciesId: TNF, TGF, IL1beta from enums.h).
   */
  float meanNeighborConcentration(SpeciesId species);

  /** Local patch concentration of a species. */
  float patchChemConcentration(SpeciesId species, int patch_index);

  /** Add to per-tick secretion accumulator for a species at patch. */
  void addPatchChemSecretion(SpeciesId species, int patch_index, float delta);

  /** Local patch chemotaxis signal value. */
  float patchChemotaxis(int patch_index);

  /** Move toward the neighbor with highest chemotaxis signal. */
  bool moveTowardChemotaxis();

  /*
   * Description:	Determines the number of cells of type cellIndex from 27
   * neighboring patches.
   *
   * Return: number of cells of type cellIndex from 27 neighboring patches.
   *
   * Parameters: cellIndex  -- enumic value of agent_t for cell type
   */
  int countNeighborCells(int cellIndex);

  /*
   * Description:	Determines the number of ECM proteins of type ECMIndex
   * from 27 neighboring patches.
   *
   * Return: number of ECM proteins of type ECMIndex from 27 neighboring
   * patches.
   *
   * Parameters: ECMIndex  -- enumic value of agent_t for ECM protein type
   */
  int countNeighborECM(int ECMIndex);

  /*
   * Description:	Move to a random neighboring patch while respecting
   * rules dictating movement between tissue types.
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
   *             dx        -- Difference in x-coordinate between current and new
   * agents. dy        -- Difference in y-coordinate between current and new
   * agents. dz        -- Difference in z-coordinate between current and new
   * agents. NOTE: dz = 0 because it is only 2D for now.
   */
  virtual void copyAndInitialize(Agent *original, int dx, int dy, int dz = 0);

  /* --------------------------------------------------------------------------
   */
  /*                              PUBLIC VARIABLES */
  /* --------------------------------------------------------------------------
   */
  int life[2]; // Number of lives remaining at the beginning and end of each
               // tick
  bool activate[2]; // Whether agent is activated or not at the beginning and
                    // end of each tick
  int color[2];     // Agent's color at the beginning and end of each tick
  float size[2];    // Agent's size at the beginning and end of each tick
  int type[2];      // Agent's type at the beginning and end of each tick

  /* --------------------------------------------------------------------------
   */
  /*                              STATIC VARIABLES */
  /* --------------------------------------------------------------------------
   */
  static BMWorld *agentWorldPtr; // Pointer from an agent to a BMWorld

  /** Chemistry facade for the current world (null before BMWorld chem init). */
  static ChemicalEnvironment *chemicalEnvironment();

  static Patch *agentPatchPtr; // Pointer from an agent to a Patch
  static ECM *agentECMPtr;     // Pointer from an agent to an ECM
  static int nx, ny, nz; // Number of patches in x,y,z dimensions of the world
  static int dX[27], dY[27],
      dZ[27]; // Difference in x-,y-,z-coordinates to neighboring patches

#ifdef MODEL_3D
  static int neighbor[27]; // Array of neighboring patches in 3D model
#else
  static int neighbor[8]; // Array of neighboring patches in 2D model
#endif

  /* --------------------------- PROTECTED VARIABLES --------------------------
   */
protected:
  int ix[2], iy[2], iz[2]; // Agent position in x,y,z dimensions at the
                           // beginning and end of each tick
  int index[2];  // Patch row major index for agent at the beginning and end of
                 // each tick
  bool alive[2]; // Life status of agent at the beginning and end of each tick
  bool realDeath[2]; // Tracking of true cell death (not false death for
                     // replacement/differentiation of cell)
  int doublings[2]; // Tracking of cell doublings (cells reach senescence, i.e.,
                    // non-proliferation, after ~100 doublings)

  abm::rng::Stream rng; // This agent's random stream

  /*
   * Description:	Opens this agent's random stream. The stream is keyed by the
   *              patch the agent is born on and the tick it is born in, which
   *              together identify it deterministically. Must be called from
   *              every agent constructor.
   *
   * Return: void
   * Parameters: birth_index -- Patch row major index the agent is born on
   */
  void openRngStream(int birth_index);
};

#endif /* AGENT_H */