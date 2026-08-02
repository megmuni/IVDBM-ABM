/*
 * parameters.h
 *
 * Calibration output helpers.
 *
 * Author: Kimberley Trickey
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 */

#ifndef PARAMETERS_H_
#define PARAMETERS_H_

#pragma once

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <iomanip>
#include <cmath>

using namespace std;
namespace util {

void outputTotalChem(BMWorld* myWorld, string filename) {
#ifdef CALIBRATION
  cout << "Outputting total chem to: " << filename << endl;
  ofstream output_file;
  output_file.open(filename.c_str(), ios::app);

  int col = 0;
  int agg = 0;
  int tp = 0;

  int nx = myWorld->nx;
  int ny = myWorld->ny;
  int nz = myWorld->nz;

  for (int in = 0; in < (nx - 1) + (ny - 1)*nx + (nz - 1)*nx*ny; in++) {
		col += myWorld->worldECM[in].ocollagen[read_t];
		col += myWorld->worldECM[in].ncollagen[read_t];
		col += myWorld->worldECM[in].fcollagen[read_t];

    agg += myWorld->worldECM[in].oaggrecan[read_t];
		agg += myWorld->worldECM[in].naggrecan[read_t];
		agg += myWorld->worldECM[in].faggrecan[read_t];
    
    tp += myWorld->worldECM[in].ocollagen[read_t];
		tp += myWorld->worldECM[in].ncollagen[read_t];
		tp += myWorld->worldECM[in].fcollagen[read_t];

		tp += myWorld->worldECM[in].oaggrecan[read_t];
		tp += myWorld->worldECM[in].naggrecan[read_t];
		tp += myWorld->worldECM[in].faggrecan[read_t];
  }

  tp += myWorld->world_total_tnf();
  tp += myWorld->world_total_tgf();
  tp += myWorld->world_total_il1beta();

  output_file << fixed << myWorld->cells.actualSize() << "\t";
  output_file << col << "\t";
  output_file << agg << "\t";
  output_file << tp << "\t";

  output_file << fixed << myWorld->world_total_tnf() << "\t";
  output_file << myWorld->world_total_tgf() << "\t";
  output_file << myWorld->world_total_il1beta() << "\t";
  
  output_file << endl;
  
#endif  // ifdef CALIBRATION
	return;
}
}  // namespace util

#endif /* PARAMETERS_H_ */
