/* 
 * Chemical.h
 * 
 * File Contents: Contains declarations for Chemical class
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *               Meghana Munipalle
 */

#ifndef Chemical_H
#define	Chemical_H

#include "../FieldVariable.h"
#include "../../common.h"
#include <vector>

using namespace std;
class World; 

/* 
 * CHEMICAL CLASS DESCRIPTION:           Chemical is a derived class of the parent class FieldVariable. 
 *                                       It contains all data members related to chemical concentration and gradients. 
 *                                       Can be used to access chemical concentrations by name.
 */
class Chemical: public FieldVariable {
 public:
    /*
     * Description:	Default Chemical constructor
     *
     * Return: void
     *
     * Parameters: void
     */
    Chemical();

    /*
     * Description:	Chemical constructor
     *
     * Return: void
     *
     * Parameters: nx  -- x-coordinate of the patch the Chemical is on
     *             ny  -- y-coordinate of the patch the Chemical is on
     *             nz  -- z-coordinate of the patch the Chemical is on
     */
    Chemical(int nx, int ny, int nz);

    /*
     * Description: Chemical destructor
     *
     * Return: void
     *
     * Parameters: void
     */
    ~Chemical();


#ifdef GPU_DIFFUSE

    /* The following variables keep track of a patch's cytokine level immediately after the GPU diffusion is performed */
    float *tTNF, *tTGF, *tIL1beta;

#endif

    // The following variables keep track of a patch's cytokine levels (pTNF, pTGF,...) and the change in a patch's cytokine levels throughout the current tick (dTNF, dTGF,...):
    float* pTNF, *dTNF;
    float* pTGF, *dTGF;
    float* pIL1beta, *dIL1beta;
    float *pcellgrad;    // Keep track of the strength of the gradients that attract cells,
    float totalTNF, totalTGF, totalIL1beta;     // Keep track of the total cytokine levels in the world
};

#endif	/* Chemical_H */