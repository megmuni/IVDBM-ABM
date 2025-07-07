/* 
 * File:   main.cpp
 * Author: alireza
 *
 * Created on May 19, 2013, 7:34 AM
 */
/*****************************************************************************
 ***  Copyright (c) 2013 by A. Najafi-Yazdi                                ***
 *** This computer program is the property of Alireza Najafi-Yazd          ***
 *** and may contain confidential trade secrets.                           ***
 *** Use, examination, copying, transfer and disclosure to others,         ***
 *** in whole or in part, are prohibited except with the express prior     ***
 *** written consent of Alireza Najafi-Yazdi.                              ***
 ******************************************************************************/

#include <cstdlib>
#include <math.h> 


using namespace std;

#include "../src/World/World.h"
#include "../src/Patch/Patch.h"
#include "../src/Agent/Agent.h"
#include "../src/ECM/ECM.h"
#include "../src/FieldVariable/FieldVariable.h"
#include "../src/enums.h"

// We do not use this anymore

void initializeFieldVariable(World* world, int ivar);

int main(int argc, char** argv) {

//    World myWorld;
    
    
    /*=================================
     Setting up the grid
     */
/*
    int nx=51;
    int ny=51;
    int nz=1;
    
    REAL x_min=-5.0;
    REAL x_max=5.0;
    REAL y_min=-10.0;
    REAL y_max=10.0;
    REAL z_min=0.0;
    REAL z_max=0.0;
   
    myWorld.setupGrid(nx,ny,nz,x_min,x_max, y_min, y_max, z_min, z_max);
    
    int ntype_fv=1; //only temperature
    myWorld.allocate_FieldVariable_Type(ntype_fv);
    
    FieldVariable temp; //defining temperature
    
    temp.BC_Type[0]=BC_Dirichlet;
    temp.BC_Type[1]=BC_Dirichlet;
    temp.BC_Type[2]=BC_Dirichlet;
    temp.BC_Type[3]=BC_Dirichlet;
    
    temp.dirichlet_bc_value[0]=1.;
    temp.dirichlet_bc_value[1]=1.;
    temp.dirichlet_bc_value[2]=1.;
    temp.dirichlet_bc_value[3]=1.;
    
    
    myWorld.field_var_type[0]=&temp; //assigning temperature 
    
    int ivar=0;
    initializeFieldVariable(&myWorld, ivar);
    
    
    const char* filename="output/init_field1.vts";
    REAL t_0=0.0;
    myWorld.outputWorld_VTK_binary(filename,t_0);
    
    
    int it_max=50;
    REAL dt=0.4;
    for (int it=0; it<it_max; it++)
    {
        myWorld.diffuse_FieldVar(dt*it,dt);
    }
    
    
    REAL t_max=dt*it_max;
    myWorld.outputWorld_VTK_binary("output/finale_field1.vts",t_max);
    return 0;
*/
}

void initializeFieldVariable(World* world, int ivar)
{
/*
     int nx = world->nx;
     int ny = world->ny;
     int nz = world->nz; 
     REAL sigma=1.5;
     REAL sigma2=sigma*sigma;
     for (int iz=0; iz<nz; iz++)
         for (int iy=0; iy<ny; iy++)
             for (int ix=0; ix<nx; ix++)
             {
                 int in=ix+iy*nx+iz*ny*nx;
                 REAL x=world->x[in];
                 REAL y=world->y[in];
                 REAL z=world->z[in];
                 
                 REAL r2=(x*x+y*y+z*z);
                 
                 
                 world->field_var[ivar][in]=1.+exp(-r2/sigma2);
             }
*/
} 

