/* 
 * File: Patch.h
 *
 * File Contents: Contains declarations for the Patch class
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */ 

#ifndef PATCH_H
#define	PATCH_H

//#include "../Agent/Agent.h" // Agent class is not needed here
#include "../FieldVariable/FieldVariable.h"

#include "../common.h"

#include <stdlib.h>
#include <vector>

using namespace std;

/*
 * PATCH CLASS DESCRIPTION:              The Patch class manages Patches (grid points in the world). 
 *                                       It is used to manage whether a patch is occupied, the patch type the patch color, 
 *                                       the type of agent on the patch, and to perform updates.
 */
class Patch {
//friend class BMWorld; 
 public:
    /*
     * Description:	Default patch constructor. Initializes patch class members.
     *
     * Return: void
     *
     * Parameters: void
     */
    Patch();

    /*
     * Description:	Patch constructor. Initializes patch class members.
     *
     * Return: void
     * 
     * Parameters: x      -- Position of patch in x dimension
     *             y      -- Position of patch in y dimension
     *             z      -- Position of patch in z dimension
     *             index  -- Patch row major index
     */
    Patch(int x, int y, int z, int index);

    /* 
     * Description:     Copy contructor
     *
     * Return: void
     */
    Patch ( const Patch &obj); 

    /*
     * Description:	Patch destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    ~Patch();

    /* 
    * Description: Overload = operator 
    */
    Patch& operator=(const Patch& originalPatch);


    /*
     * Description:	Determines whether the patch is occupied at the beginning of the tick. 
     *              Is used for atomic-move support.
     *
     * Return: True is this patch is occupied	(According to the read-entry)
     *
     * Parameters: void
     */
    bool isOccupied();

    /*
     * Description:	Determines whether the patch is occupied at the end of the tick.
     *
     * Return: True is this patch is occupied (According to the write-entry)
     *
     * Parameters: void
     */
    bool isOccupiedWrite();

    /*
     * Description:	Test and set the write-entry of the occupancy of this patch 
     *              and set the dirty flag so an update function gets called for this patch.
     *
     * Return: True if patch was already set (According to the write-entry)
     *
     * Parameters: void
     */
    bool setOccupied();

    /*
     * Description:	Test and set the read-entry of the occupancy of this patch
     *
     * Return: True if patch was already set (According to the read-entry)
     *
     * Parameters: void
     */
    bool setOccupiedRead();

    /*
     * Description:	Clear the write-entry of the occupancy of this patch and 
     *              set the dirty flag so an update function gets called for this patch.
     *
     * Return: void
     *
     * Parameters: void
     */
    void clearOccupied();

    /*
     * Description:	Clear the read-entry of the occupancy of this patch
     *
     * Return: void
     *
     * Parameters: void
     */
    void clearOccupiedRead();

    /*
     * Description:	Test and set the write-entry of the occupancy of this patch
     * 				and set the dirty flag so an update function gets called for this patch
     *
     * Return: void
     *
     * Parameters: void
     */
    void setOccupiedLight();

    /*
     * Description:	Test and set the read-entry of the occupancy of this patch
     * 				and set the dirty flag so an update function gets called for this patch
     *
     * Return: void
     *
     * Parameters: void
     */
    void setOccupiedReadLight();

    /*
     * Description:	Outputs patch coordinates and type (damage, unidentifiable)
     *
     * Return: void
     *
     * Parameters: void
     */
    void lookUpType();				

    /*
     * Description:	Determines the color of the patch.
     *
     * Return: Enumic value of color_t for patch color
     *
     * Parameters: void
     */
    int getColorfromType(); 

    /*
     * Description:	Determines whether there is an agent of type 'agentType' on the current patch. 
     *
     * Return: True if there is agent of type 'agentType' on current patch
     *
     * Parameters: agentType  -- integer definining the type of agent 
     *                           (0: platelet, 1: chondrocytes, 2: macrophages,
     *                            3: neutrophils, 4: original collagen,
     *                            5: new collagen, 6: fragmented collagen
     *                            7: original aggrecan, 8: new aggrecan, 
     *                            9: fragmented aggrecan, 
     *                            10: original hyaluronan, 11: new hyaluronan,
     *                            12: fragmented hyaluronan)
     */
    bool agentOnPatch(int agentType); 

    /*
     * Description:	Update the current patch with all changes in occupation, health, and damage.
     *
     * Return: void
     *
     * Parameters: void
     */
    void updatePatch();

    /*
     * Description:	Update the color attribute of the current patch to an agent color if applicable, or its tissue type color. 
     *
     * Return: void
     *
     * Parameters: void
     */
    void render();

    /*************************************************************************
     * CONSTANT ATTRIBUTES                                                   *
     *************************************************************************/
    int indice[3];    // Used to store the patch's location in x,y,z dimensions of world
    int index;        // Used to store the patch row major index
    bool inDamzone;   // Whether the patch is part of the damaged zone of tissue
    bool initHA;      // Whether the patch can be a center for sprouting original hyaluronan
    bool initcollagen;   // Whether the patch can be a center for sprouting original collagen
    bool initaggrecan;    // Whether the patch can be a center for sprouting original aggrecan

    /*************************************************************************
     * VARIABLE ATTRIBUTES                                                   *
     *************************************************************************/
    int type[2];     // Keeps track of type of tissue on patch at beginning and end of each tick
    int color[2];    // Keeps track of the color of the patch (either from an agent or tissue)  at the beginning and end of each tick
    float health[2];           // Keeps track of the health of the patch (0 or 100) at the beginning and  end of each tick 
    unsigned int damage[2];    // Keeps track of the damage on the patch (0 to the max number of fragmented ECM proteins) at the beginning and end of each tick
    int occupiedby[3];    // Keeps track of the type of agent in the patch (or holds 'unoccupied') at the beginning and end of each tick
    bool dirty;           // Whether the patch has been modified and needs updating

    /*************************************************************************
     * STATIC ATTRIBUTES                                                     *
     *************************************************************************/
    /* Keeps track of the number of patches of each type. 
     *      Indexed according to the enumic values of patches_t 
     *      (0: blood, 1:tissue, 2:epithelium, 3:capillary, 4:damage, 5:CaAlg)
     */
    #ifdef MODEL_SCAFFOLD
        static float numOfEachTypes[6];
    #else
        static float numOfEachTypes[5];
    #endif
    
    /*************************************************************************
     * The following four functions were part of the effort to "standardize" *
     * interaction of patches with agents. Consider finishing these if want  *
     * to implement more than 1 agent per patch. (probe, moveOnPatch, bornOn *
     * moveOffPatch, are already implemented in Agent functions, eg, move.)  *
     * int probe(int cellType);                                              *
     * void moveOnPatch(Agent* incomingCell);                                *
     * void moveOffPatch(Agent* leavingCell);                                *
     * void bornOn(Agent* newBorn);                                          *
     * vector<Agent*> patchAgentPtr;                                         *
     *************************************************************************/

 private:
    /* This atomic_flag variable is used to ensure thread-safety.
     * Since there's a chance two threads might try to move their cells to the
     * same patch at the same time, all patch availability updates need to be
     * done atomically.
     *
     * The boolean variable is here for a light read since std::atomic_flag 
     * doesn't support a read without writing, this this boolean variable would
     * help reduce the work when we just want to perform a read.
     *
     * NOTE: The atomic_flag and the boolean variables need to be synchronized
     * at all times. Thus, it's crucial they are kept private and ALL writes to
     * patch occupancies need to be done via a setter function. */
    atomic_flag* occupied[2];
    // Whether the patch is occupied at the beginning and end of each tick
    bool occupied_bool[2];
};

#endif	/* PATCH_H */