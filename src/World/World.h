/*
 * File: World.h 
 *
 * File Contents: Contains declarations for the World class
 *   
 * Author: alireza
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#ifndef WORLD_H
#define	WORLD_H

#include "../Agent/Agent.h"
#include "../ArrayChain/ArrayChain.h"
#include "../FieldVariable/FieldVariable.h"
#include "../Patch/Patch.h"

#include <stdlib.h>
#include <vector>

using namespace std;
//extern class dev_Patch;

/*
 * WORLD CLASS DESCRIPTION:          The World class manages all worlds in the model. 
 *                                   It is used to set up the dimensions of the world.
 */
class World {
 public:
    /*
     * Description:	Default World constructor. 
     *
     * Return: void
     *
     * Parameters: void
     */
    World();

    /*
     * Description:	World constructor. 
     *
     * Return: void
     *
     * Parameters: orig  -- Pointer to an original World
     */
    World(const World& orig);

    /*
     * Description:	Virtual World destructor.
     *
     * Return: void
     *
     * Parameters: void
     */
    virtual ~World();

    /*
     * Description:	Initializes the grid size and dimensions
     *
     * Return: void
     *
     * Parameters: nx, ny, nz    -- Number of grid points (patches) in x,y,z dimensions
     *             x_min, x_max  -- Min and max coordinates in x
     *             y_min, y_max  -- Min and max coordinates in y
     *             z_min, z_max  -- Min and max coordinates in z
     */
    void setupGrid(int nx, int ny,int nz, REAL x_min, REAL x_max, REAL y_min, REAL y_max, REAL z_min, REAL z_max); //!< set up the dimensions

    /*
     * Description:	Write vtk output file "filename" for animation
     *
     * Return: void
     *
     * Parameters: filename  -- Path to output file
     *             t         -- time
     */
    void outputWorld_VTK_binary(const char* filename, REAL t);
        int nx, ny, nz;     // Number of grid points (lattices) (patches) on this world in x,y,z
        REAL dx, dy, dz;    // Mesh (patch) size in x,y,z dimenstions
        REAL x_min, x_max;    // Max and min coordinate in x 
        REAL y_min, y_max;    // Max and min coordinate in y 
        REAL z_min, z_max;    // Max and min coordinate in z 
        vector<REAL> x,y,z;    // Spatial coordinates of the grid
        vector<vector<REAL> > field_var;      // Two dimensional array of field variables 
        unsigned seed;    // For generating random numbers
};

#endif	/* WORLD_H */