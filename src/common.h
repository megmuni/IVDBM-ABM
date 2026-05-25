/*
 * File: common.h
 *
 * File Contents: Preprocessor directives for adjusting how the model executes.
 *
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#ifndef COMMON_H_
#define COMMON_H_
#ifdef __GNUC__
	#include <atomic>
	//#include <cstdatomic>
#else
	#include "atomic"
#endif
#include "enums.h"

/*****************************************************************************
 * SENSITIVITY ANALYSIS AND  CALIBRATION                                     *
 *****************************************************************************/
// Modify normal values of parameters for S.A. or calibration purposes
#define CALIBRATION 	

// Used to print new values of parameters during S.A. or calibration
#define PRINT_PARAMETER_VALUES


/*****************************************************************************
 * PROFILING FLAGS - FOR TIMING INFORMATION                                  *
 *****************************************************************************/
// Profile major sections of go():
//   //#define PROFILE_MAJOR_STEPS  // Turned this off for SA

// Profile each cell type:
//  //#define PROFILE_CELL_FUNC  // Turned this off for SA

// Profile each section of chemical diffusion:
#define PROFILE_CHEM_DIFF_DETAILED

// Profile each section of cell seeding:
#define PROFILE_CELL_PROLIF_DETAILED

// Profile time spent on chemical diffusion by each thread:
//#define PROFILE_THREAD_LEVEL_CHEM_DIFF

// Profile each section of ECM updates:
//    //#define PROFILE_ECM_UPDATE  // Turned this off for SA

// Profile each section of patch updates
//#define PROFILE_PATCH_UPDATE  // Turned this off for SA

/*****************************************************************************
 * PARALLEL OPTION FLAGS                                                     *
 *****************************************************************************/
//#define PAR_BY_CELLTYPE // Parallel by cell types	(Obsolete)

/*****************************************************************************
 * OPTIMIZATION FLAGS                                                        *
 *****************************************************************************/
//#define VECTORIZE // TODO(Nuttiiya)
//#define OPT_CELL_SEEDING // TODO(Nuttiiya)
//#define ECM_UNROLL_LOOP // TODO(Nuttiiya)

/* With OPT_ECM defined, we assume that we access ECM::HAlife in this order:
 * 1. Decrement life by calling decrement(int n) on HAlife in ECMFunction()
 * 2. Determine number of HA by calling HAlife.size() in fragmentHA()
 * 3. Remove dead HAs in HAlife in updateECM()
 * NOTE: Steps 1 & 2 can occur in either order but must preceed step 3. */
#define OPT_ECM

/*****************************************************************************
 * OPTION FLAGS                                                              *
 *****************************************************************************/
#define BIOMARKER_OUTPUT // Enable Biomarker Output
//#define PARAVIEW_RENDERING // Enable Paraview rendering of Patches and Cells
//#define PARAVIEW_ECM_CHEM // Enable Paraview rendering of ECM and Chem
#define MODEL_3D // Run in 3D

// Initialize Alg/Ca Scaffold (for in vitro case define MODEL_SCAFFOLD, in vivo case define MODEL_VOCALFOLD && MODEL_SCAFFOLD)
#define MODEL_SCAFFOLD
//#define PEVOC_SCALE // Use PEVOC scale parameters 	(Obsolete?)
#define PDE_DIFFUSE // Use PDE based chemical diffusion
//#define COLLECT_CELL_INS_DEL_STATS // TODO(Nuttiiya)

#define NUM_TICKS	240
#define NUM_THREAD 16 		// Number of threads for parallelization
#define MAX_NUM_THREADS 16  // Maximum number of threads for parallelization
#define DEFAULT_TID	0 		// Default thread identification number
#define DEFAULT_DATA_SMALL 1<<15  // Constant for ArrayChain. Array size: ~ 32k (32768)
#define DEFAULT_DATA_MEDIUM 1<<18 // Constant for ArrayChain. Array size: ~262k (262144)
#define DEFAULT_DATA_LARGE 1<<19  // Constant for ArrayChain. Array size: ~524k (524288)

/* NOTE: Please ignore these options (V_a, V_b) since they are for the currently incorrect diffuseChem() */
//#define V_a
#define V_b
//#define PRINT_KERNEL

// TODO(Nuttiiya)
#ifdef VECTORIZE
typedef float v4sf __attribute__ ((vector_size(sizeof(float)*4)));
union f4vector
{
	v4sf v;
	float f[4];
};
#endif

typedef int sizeType; // TODO(Nuttiiya)

/* Description:	Check if a field is modified in this tick (i.e. read and write entry are not the same). Used for optimizations in updates
 *
 * Return: True if modified
 * Parameters: arr  -- Pointer to attribute to check
 */
template <typename T>
bool isModified(T* arr){
	return arr[read_t] != arr[write_t];
}
#endif /* COMMON_H_ */