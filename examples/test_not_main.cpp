/* 
 * File: test_not_main.cpp
 *
 * File Contents: Contains main method of model and various output functions.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

// Include C/C++ libraries
#include <cstdlib>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <fstream>
#include <cstring>
#include <sstream>
#include <time.h>
#include <sys/time.h>

//Include local libraries
#include "../src/ArrayChain/ArrayChain.h"
#include "../src/World/World.h"
#include "../src/World/Usr_World/biomaterialWorld.h"
#include "../src/Agent/Agent.h"
#include "../src/Agent/Usr_Agents/Cell.h"
#include "../src/FieldVariable/FieldVariable.h"
#include "../src/Patch/Patch.h"
#include "../src/ECM/ECM.h"
#include "../src/enums.h"
#include "../src/FieldVariable/Usr_FieldVariables/Chemical.h"
#include "../src/Utilities/input_utils.h"
#include "../src/Utilities/parameters.h"
#include "../src/World/world_vtk_export.h"

using namespace std;

/*
 * Description:	Outputs the color of each patch to the given file.
 *              Outputs the color of the agent that is one the patch and if there is no agent, then outputs the color of the patch type.
 *              Also prints the number of patches that are chondrocyte, new chondrocyte, black, and the total number of patches.
 *
 * Return: 0 if succesful
 *
 * Parameters: BMWorld*  -- Pointer to the wound healing world whose patches' colors should be outputted.
 *             char*     -- Output file name
 */
int outputColor(BMWorld*,char*);

/*
 * Description:	Outputs the color of each patch to the given file.
 *              Outputs the color of the ECM that is one the patch and if there is no ECM, then outputs the color of the patch type.
 *
 * Return: 0 if succesful
 *
 * Parameters: myWorld   -- Pointer to the wound healing world whose patches' colors should be outputted.
 *             fileName  -- Output file name
 */
int outputECM(BMWorld* myWorld, char* fileName);

/*
 * Description:	Outputs the chemical concentration of the given chemical on each patch to the given file.
 *
 * Return: 0 if succesful
 *
 * Parameters: BMWorld*  -- Pointer to the wound healing world whose patches' colors should be outputted.
 *             char*     -- Output file name
 *             int       -- Enumic value of the chemical type to output
 */
int outputChem(BMWorld*, char*, int);

/*
 * Description:	Main method for model simulation. Sets up the world and executes each tick of the simulation. 
 * 				Calculates cell insertion and deletion statistics as well as execution time information.
 *
 * Return: 0 if succesful
 *
 * Parameters: argc  -- Number of command-line arguments
 *             argv  -- Command-line arguments
 */
int main(int argc, char** argv) {

/* -------------------------------------------------------------------------- */
/*                        SET UP A WOUND HEALING WORLD                        */
/* -------------------------------------------------------------------------- */
	// Get baseline cell and chemical values that user specified:
	util::processOptions(argc, argv);
	util::ensureOutputDir();
	util::writeRunParamsJson(argc, argv);
	util::printOptions();
	util::processParameters("Sample.txt");
	clock_t tStart = clock();
	BMWorld myWorld(util::getWorldXWidth(), util::getWorldYWidth(), util::getWorldZWidth(), util::getPatchWidth());
	myWorld.outputWorld_csv();
	//printf("Setup Execution time: %.2fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);

/* -------------------------------------------------------------------------- */
/*              EXECUTE EACH TICK (30 MINUTES EACH) OF THE MODEL              */
/* -------------------------------------------------------------------------- */
	int numTicks = util::getNumTicks();
	struct timeval start, end; 		// Timing information
	long elapsed_times[numTicks];   // in milliseconds

	for (int tick = 0; tick < numTicks; tick++) {
		//myWorld.debugInfo(); // debug function
		if (util::paraviewEnabled() && tick % 12 == 0) {
			char paraview_dir[512];
			util::ensureOutputSubpath("paraview");
			util::makeOutputPath(paraview_dir, sizeof(paraview_dir), "paraview");
			export_world_timestep(myWorld, tick, paraview_dir);
		}

		// Run the simulation for 1 tick (30 min):
		clock_t t1 = clock();
		cerr << "executing go() at tick " << tick << " ..." << endl;
		cout << "entering go() at tick " << tick << endl;
		gettimeofday(&start, NULL);
		myWorld.go();

		#ifdef BIOMARKER_OUTPUT
			myWorld.outputWorld_csv();
		#endif
		
		gettimeofday(&end, NULL);
		elapsed_times[tick] = (end.tv_sec*1000 + end.tv_usec/1000) - (start.tv_sec*1000 + start.tv_usec/1000); // Store execution time for this tick
		//cout << "	this go() execution took " << elapsed_times[tick] << " ms" << endl;
	}

	#ifdef CALIBRATION // Used for sensitivity analysis
		char calib_path[512];
		util::ensureOutputSubpath("SensitivityAnalysis");
		util::makeOutputPath(calib_path, sizeof(calib_path),
		                     "SensitivityAnalysis/FinalTotalChemVR.dat");
		util::outputTotalChem(&myWorld, calib_path);
	#endif

#ifdef COLLECT_CELL_INS_DEL_STATS
	int addedSum = 0;
	int erasedSum = 0;
	int addedMax = myWorld.numAddedCells[0];
	int erasedMax = myWorld.numErasedCells[0];
	int addedMin = myWorld.numAddedCells[0];
	int erasedMin = myWorld.numErasedCells[0];
	int addedMin_NoZero = myWorld.numAddedCells[0];
	int erasedMin_NoZero = myWorld.numErasedCells[0];

 /* -------------------------------------------------------------------------- */
 /*                          CELL INSERTION STATISTICS                         */
 /* -------------------------------------------------------------------------- */
	cout << "	Cell Insertion Stats: " << endl;
	for (int i = 0; i < numTicks; i++) {
	  	cout << "		tick " << i << ": "<< endl;
		cout << "		# Added Chonds:	   "<< myWorld.numAddedChondros[i] << endl;
		cout << myWorld.numAddedCells[i] << endl;
		addedSum += myWorld.numAddedCells[i];
		addedMax = myWorld.numAddedCells[i] > addedMax? myWorld.numAddedCells[i] : addedMax;
		addedMin = myWorld.numAddedCells[i] < addedMin? myWorld.numAddedCells[i] : addedMin;
		addedMin_NoZero = (myWorld.numAddedCells[i] < addedMin) && (myWorld.numAddedCells[i] > 0)? myWorld.numAddedCells[i] : addedMin;
	}

	cout << "		Min:		" << addedMin << endl;
	cout << "		Min (no 0):	" << erasedMin << endl;
	cout << "		Max:		" << addedMax << endl;
	cout << "		Average:	" << addedSum/numTicks << endl;
	cout << endl;

 /* -------------------------------------------------------------------------- */
 /*                          CELL DELETION STATISTICS                          */
 /* -------------------------------------------------------------------------- */
	cout << "	Cell Deletion Stats: " << endl;
	for (int i = 0; i < numTicks; i++) {
		cout << "		tick " << i << ": "<< endl;
		cout << "	    # Erased Chonds:	" << myWorld.numErasedChondros[i] << endl;
		cout << myWorld.numErasedCells[i] << endl;
		erasedSum += myWorld.numErasedCells[i];
		erasedMax = myWorld.numErasedCells[i] > addedMax? myWorld.numErasedCells[i] : addedMax;
		erasedMin = myWorld.numErasedCells[i] < addedMin? myWorld.numErasedCells[i] : addedMin;
		erasedMin_NoZero = (myWorld.numErasedCells[i] < addedMin) && (myWorld.numErasedCells[i] > 0)? myWorld.numErasedCells[i] : addedMin;
	}

	cout << "		Min:		" << erasedMin << endl;
	cout << "		Min (no 0):	" << erasedMin << endl;
	cout << "		Max:		" << erasedMax << endl;
	cout << "		Average:	" << erasedSum/numTicks << endl;
#endif

/* -------------------------------------------------------------------------- */
/*                             TIMING INFORMATION                             */
/* -------------------------------------------------------------------------- */
	int acc = 0;
	for (int i = 0; i < numTicks; i++) {
		//cout << "Tick " << i << " took " << elapsed_times[i] << " ms" << endl;
		acc += elapsed_times[i];
	}
	//cout << "Average time:	" << acc/numTicks << " ms" << endl;

	util::writeRunTimingJson(
	    acc, numTicks, (double)(clock() - tStart) / CLOCKS_PER_SEC);

	return 0;
} // End Main

int outputColor(BMWorld* myWorld, char* fileName) {
	int dam = 0, tissue = 0, numCaAlg = 0, cells = 0, newCell = 0, stems = 0, progens = 0, nps = 0, total = 0, black = 0, actDam = 0;
	ofstream outfile(fileName);

  /* Prepare legacy VTK file format for for visualization with Paraview 3.0 */
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Really cool data " << endl;
	outfile << "ASCII " << endl;
	outfile << "DATASET STRUCTURED_POINTS " <<endl;
	outfile << "DIMENSIONS " << myWorld->nx << " " << myWorld->ny << " " << myWorld->nz << endl;
	outfile << "ORIGIN 0 0 0 " << endl;
	outfile << "SPACING 1 1 1 " << endl;
	outfile << "POINT_DATA " << myWorld->nx*myWorld->ny*myWorld->nz << endl;
	outfile << "SCALARS Color int 1 " << endl;
	outfile << "LOOKUP_TABLE default " << endl;

	int in = 0;
	for (int iz = 0; iz < myWorld->nz; iz++) {
		for (int iy = 0; iy < myWorld->ny; iy++) {
			for (int ix = 0; ix < myWorld->nx; ix++) {

				if (ix == (myWorld->nx - 1) && iy == (myWorld->ny - 1)) {
					outfile << "195";  // Visualization color legend upper bound

				} else if (ix == (myWorld->nx - 2) && iy == (myWorld->ny - 1)) {
					outfile << "0 ";  // Visualization color legend lower bound

				} else {
					in = ix + iy*myWorld->nx + iz*myWorld->nx*myWorld->ny;
					outfile << (myWorld->worldPatch[in].color[read_t]) << " "; // Output the color on each patch

					// Count the number of cells and patch types:
					if (myWorld->worldPatch[in].occupiedby[read_t] == cell)
						cells++;

					else if (myWorld->worldPatch[in].type[read_t] == stem)
						stems++;

					else if (myWorld->worldPatch[in].type[read_t] == progen)
						progens++;

					else if (myWorld->worldPatch[in].type[read_t] == np)
						nps++;

					else if (myWorld->worldPatch[in].type[read_t] == CaAlg) 
						numCaAlg++;

					else if (myWorld->worldPatch[in].type[read_t] == damage)
						black++;
					else {
						//cout << "the color is " << myWorld->worldPatch[in].color << " and location is ";
						//cout << ix << " " << iy << endl;
					}
					total++;
				}
			}
			outfile << endl;
		}
	}

    /* // Output the cell & patch type counts: 
    cout << "file name is " << fileName << endl;
	cout << " the counts are: all cells " << cells << " , stems " << stem << " , progenitors " << progens << " , NP cells " << nps << endl;		
	cout << " damage " << black << " , CaAlg " << numCaAlg << " , total " << total << endl;
	cout << " to check, the num of cells are " << Cell::numOfCells << endl; 
	*/
	return 0;

} // End outputColor

int outputECM(BMWorld* myWorld, char* fileName) {
  /* Prepare legacy VTK file format for for visualization with Paraview 3.0 */
	ofstream outfile(fileName);
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Really cool data " << endl;
	outfile << "ASCII " << endl;
	outfile << "DATASET STRUCTURED_POINTS " << endl;
	outfile << "DIMENSIONS " << myWorld->nx << " " << myWorld->ny << " " << myWorld->nz << endl;
	outfile << "ORIGIN 0 0 0 "<< endl;
	outfile << "SPACING 1 1 1 "<< endl;
	outfile << "POINT_DATA " << myWorld->nx*myWorld->ny*myWorld->nz << endl;
	outfile << "SCALARS Color int 1 " << endl;
	outfile << "LOOKUP_TABLE default " << endl;

  // Assign and output the appropriate color to each patch:
    int in = 0; 
	for (int iz = 0; iz < myWorld->nz; iz++) {
		for (int iy = 0; iy < myWorld->ny; iy++) {
			for (int ix = 0; ix < myWorld->nx; ix++) {
                in = ix + iy*myWorld->nx + iz*myWorld->nx*myWorld->ny;

                if (myWorld->worldECM[in].empty[read_t] == false) {
					if (myWorld->worldECM[in].oaggrecan[read_t] != 0 || myWorld->worldECM[in].naggrecan[read_t] != 0) {
						outfile << caggrecan << " ";

					} else if ( myWorld->worldECM[in].faggrecan[read_t] != 0) {
						outfile << cfaggrecan << " ";

					} else if (myWorld->worldECM[in].ocollagen[read_t] != 0 || myWorld->worldECM[in].ncollagen[read_t] != 0) {
						outfile << ccollagen << " ";

					} else if ( myWorld->worldECM[in].fcollagen[read_t] != 0) {
						outfile << cfcollagen << " ";
					} 
				} else outfile << myWorld->worldPatch[in].getColorfromType() << " ";
			}
			outfile << endl;
		}
	}
	return 0;

} // End outputECM

int outputChem(BMWorld* myWorld, char* fileName, int chemIndex) {

  /* Prepare legacy VTK file format for for visualization with Paraview 3.0 */
    ofstream outfile(fileName);
	outfile << "# vtk DataFile Version 2.0" << endl;
	outfile << "Really cool data " << endl;
	outfile << "ASCII " << endl;
	outfile << "DATASET STRUCTURED_POINTS " << endl;
	outfile << "DIMENSIONS " << myWorld->nx << " " << myWorld->ny << " " << myWorld->nz << endl;
	outfile << "ORIGIN 0 0 0 " << endl;
	outfile << "SPACING 1 1 1 " << endl;
	outfile << "POINT_DATA " << myWorld->nx*myWorld->ny*myWorld->nz << endl;
	outfile << "SCALARS Color float 1 " << endl;
	outfile << "LOOKUP_TABLE default " << endl;

  // Output the chemical concentration on each patch:
	for (int iz = 0; iz < myWorld->nz; iz++) {
		for (int iy = 0; iy < myWorld->ny; iy++) {
			for (int ix = 0; ix < myWorld->nx; ix++) {
				int in = ix + iy*myWorld->nx + iz*myWorld->nx*myWorld->ny;
				const float *grid =
				    myWorld->chemical_environment()->channel_grid(chemIndex);
				outfile << grid[in] << " ";
			}
			outfile << endl;
		}
	}
	return 0;

} // End outputChem