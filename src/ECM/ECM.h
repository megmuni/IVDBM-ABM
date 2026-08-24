/*
 * File: ECM.h
 *
 * File Contents: Contains declarations for ECM class
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#ifndef ECM_H
#define	ECM_H

#include <stdlib.h>
#include <vector>
#include "../Agent/Agent.h"
#include "../FieldVariable/FieldVariable.h"
#include "../enums.h"
#include "../Utilities/rng.h"

class World;
class BMWorld;
class Patch;

/*
 * ECM CLASS DESCRIPTION:            The ECM class manages ECM proteins. 
 *                                   It is used to fragment ECM proteins, and signal and repair damage. 
 */
class ECM {
 public:
    /*
     * Description:	Default ECM constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    ECM();

    /*
     * Description:	ECM constructor. Initializes ECM class members.
     *
     * Return: void
     * 
     * Parameters: x      -- Position of ECM manager in x dimension
     *             y      -- Position of ECM manager in y dimension
     *             z      -- Position of ECM manager in z dimension
     *             index  -- Patch row major index of ECM manager in world
     */
    ECM(int x, int y, int z, int index); 

    /*
     * Description:	ECM destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    ~ECM();

    /*
     * Description:	Performs biological function of ECM proteins.
     *
     * Return: void
     *
     * Parameters: void
     */
    void ECMFunction();

    /*
     * Description:	Repairs the damage on the patch where the ECM manager is.
     *
     * Return: void
     *
     * Parameters: void
     */
    void repairDamage();

    /*
     * Description:	Creates damage on the patch where the ECM manager is. 
     *              If there are activated macrophages or activated neutrophils already on the patch, the damage is cleared.
     *
     * Return: void
     *
     * Parameters: void
     */
    void dangerSignal();

    /*
     * Description:	Clears all flags present in the current tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void resetrequests();

    /*
     * Description:	Breaks down each new collagen ECM protein on this ECM 
     *              manager's patch into two fragmented collagen ECM proteins
     *
     * Return: void
     *
     * Parameters: void
     */
    void fragmentNCollagen();

    /*
     * Description:	Breaks down each new aggrecan ECM protein on this ECM 
     *              manager's patch into two fragmented aggrecan ECM proteins
     *
     * Return: void
     *
     * Parameters: void
     */
    void fragmentNAggrecan();

    /*
     * Description:	Breaks down each hyaluronan ECM protein on this ECM 
     *              manager's patch into two fragmented hyaluronan ECM proteins
     *
     * Return: void
     *
     * Parameters: void
     */
    void fragmentHA();

    /*
     * Description:	Update given ECM manager with all requests for fragment ECM from neighboring patches
     *
     * Return: void
     *
     * Parameters: void
     */
    void updateECM();

    /* Description: Assigns true to the ECM manager's empty[] array at the end 
     *              of the tick if there are no ECM proteins for this manager.
     * 
     * Return: void
     * 
     * Parameters: void
     */
    void isEmpty();

#ifdef OPT_ECM
    /*
     * Description:	Set dirty flag to indicate change in this ECM's attributes.
     *              Should be called whenever this ECM manager modifies its own attributes. 
     *              Dirty flag is checked when ECMupdate() is called to see if update is necessary.
     *
     * Return: void
     *
     * Parameters: void
     */
    void set_dirty();

    /*
     * Description:	Set request-dirty flag to indicate that this ECM manager has requested ECM fragments on a neighbor.
     *              Should be called whenever this ECM manager submits a request for fragments on a neighboring patch.
     * 				Request-dirty flag is checked when reset_request() is called to see if memset() on requests is necessary.
     *
     * Return: void
     *
     * Parameters: void
     */
    void set_request_dirty();

    /*
     * Description:	Set dirty-from-neighbor flag to indicate that a neighbor has requested ECM fragments from this ECM manager.
     *              Should be called whenever this ECM manager submits a request for fragments to a neighboring patch as follows:
     * 					
     *                 this->ECMWorldPtr->worldECM[neighbor_index].set_dirty_from_neighbors();
     *
     *              Dirty-from-neighbor flag is checked when updateECM() is called to see if there's a fragment request to be processed
     *
     * Return: void
     *
     * Parameters: void
     */
    void set_dirty_from_neighbors();

    /*
     * Description:	Reset dirty flag. Should be called after done using this flag in this tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void reset_dirty();

    /*
     * Description:	Reset request-dirty flag. Should be called after done using this flag in this tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void reset_request_dirty();

    /*
     * Description:	Reset dirty-from-neighbor flag. Should be called after done using this flag in this tick.
     *
     * Return: void
     *
     * Parameters: void
     */
    void reset_dirty_from_neighbors();
#endif

    /*
     * Description:	Decreases the input parameter by 1. Used to decrease the life of ECM hyaluronan proteins.
     *
     * Return: void
     *
     * Parameters: n  -- An integer for the number of lives for hyaluronan.
     */
    static void decrement(int n);

    /*************************************************************************
     * CONSTANT VARIABLES                                                    *
     *************************************************************************/
    int indice[3]; // Used to store the ECM's location in x,y,z dimensions of world
    int index;     // Used to store the patch row major index of the ECM manager.

    /*************************************************************************
     * NON-CONSTANT VARIABLES                                                *
     *************************************************************************/
    bool empty[2];     // Whether there is ECM or not at the beginning and end of each tick
    int ocollagen[2], ncollagen[2], fcollagen[2]; // The number of original collagen (tropocollagen monomer), new collagen(collagen) and fragmented collagen at the beginning and end of each tick
    int oaggrecan[2], naggrecan[2], faggrecan[2];    // The number of original aggrecan (monomer), new aggrecan and fragmented aggrecan at the beginning and end of each tick
    int HA[2], fHA[2];                            // The number of hyaluronan and fragmented hyaluronan at the beginning and end of each tick 

    int requestfcollagen[27];   // Keeps track of all neighbors' requests for fragmented collagen.
    int requestfaggrecan[27];    // Keeps track of all neighbors' requests for fragmented aggrecan.
    int requestfHA[27];         // Keeps track of all neighbors' requests for fragmented hyaluronan.

#ifdef OPT_ECM
    /* Number of lives for each hyaluronan (HA) remaining in this tick.
     * With OPT_ECM defined, we assume that we access HAlife in this order:
     * 1. Decrement life by calling decrement(int n) on HAlife in ECMFunction()
     * 2. Determine number of HA by calling HAlife.size() in fragmentHA()
     * 3. Remove dead HAs in HAlife in updateECM()
     * NOTE: Steps 1 & 2 can occur in either order but must preceed step 3.
     */
    std::vector<int> HAlife;
#else
    std::vector<int> HAlife[2];    // Number of lives for hyaluronan at the beginning and end of each tick.
#endif

    bool fcollDangerSignal[2];    // Whether or not there is a fragmented collagen signalling danger at the beginning and end of each tick 
    bool faggDangerSignal[2];    // Whether or not there is a fragmented aggrecan signalling danger at the beginning and end of each tick
    bool fHADangerSignal[2];   // Whether or not there is a fragmented hyaluronan signalling danger at the beginning and end of each tick 
    bool scarIndex[2];    // Whether there is enough collagen to form a scar at the beginning and end of each tick

#ifdef OPT_ECM
    bool dirty;                 // Whether there has been a change in this ECM's attributes    
    bool request_dirty;         // Whether this ECM manager has requested ECM fragments on a neighbor
    bool dirty_from_neighbors;  // Whether a neighbor has requested ECM fragments from this ECM manager

#endif

    /*************************************************************************
     * STATIC VARIABLES                                                      *
     *************************************************************************/
    static Patch* ECMPatchPtr;     // Pointer to array of all patches in the world
    static BMWorld* ECMWorldPtr;   // Pointer to the wound healing world
    static int maxcollagen;    // The maximum amount of collagen (of all types) allowed
    static int maxaggrecan;     // The maximum amount of aggrecan (of all types) allowed
    static int maxHA;     // The maximum amount of hyaluronan (of all types) allowed
    static int dx[27];    // x dimension displacement to each neighbor
    static int dy[27];    // y dimension displacement to each neighbor
    static int dz[27];    // z dimension displacement to each neighbor
    static int d[27];     // Index of each neighbor in dx, dy, dz arrays

    abm::rng::Stream rng;  // This ECM manager's random stream
};

#endif	/* ECM_H */