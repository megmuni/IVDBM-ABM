/* 
 * WHChemical.h
 * 
 * File Contents: Contains declarations for WHChemical class
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#ifndef WHChemical_H
#define	WHChemical_H

#include "../FieldVariable.h"
#include "../../common.h"
#include <vector>

using namespace std;
class World; 

/* 
 * WHCHEMICAL CLASS DESCRIPTION:         WHChemical is a derived class of the parent class FieldVariable. 
 *                                       It contains all data members related to chemical concentration and gradients. 
 *                                       Can be used to access chemical concentrations by name.
 */
class WHChemical: public FieldVariable {
 public:
    /*
     * Description:	Default WHChemical constructor
     *
     * Return: void
     *
     * Parameters: void
     */
    WHChemical();

    /*
     * Description:	WHChemical constructor
     *
     * Return: void
     *
     * Parameters: nx  -- x-coordinate of the patch the WHChemical is on
     *             ny  -- y-coordinate of the patch the WHChemical is on
     *             nz  -- z-coordinate of the patch the WHChemical is on
     */
    WHChemical(int nx, int ny, int nz);

    /*
     * Description: WHChemical destructor
     *
     * Return: void
     *
     * Parameters: void
     */
    ~WHChemical();


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

#endif	/* WHChemical_H */