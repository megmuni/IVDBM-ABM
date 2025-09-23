/*
 * parameters.h
 *
 * File Contents: TODO
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

/*
 * Description:	Initializes parameters from the given input file.
 *
 * Return: void
 *
 * Parameters: filename  -- Name of the file containing the parameter values
 */
void processParameters(string filename) {
#ifdef CALIBRATION
  cout << "Processing Parameters..." << endl;
  //ifstream input_file("Sample.txt", ios::in);
  //ifstream input_file;
  //input_file.open(filename.c_str(), ios::in);
   // ifstream input_file(filename, ios::in);
  //  ifstream input_file("NoExists.txt", ios::in);


    ifstream input_file(filename);

    if (!input_file) {
        cerr << "Could not open input file for processing parameters" << endl;
        return;
    }

    string line;
    if (!getline(input_file, line)) {
        cerr << "Failed to read data line from file" << endl;
        return;
    }

    stringstream lineStream(line);
    string value;

    /* // Read and print each value for debugging
      int index = 0;
      while (getline(lineStream, value, '\t')) {
          cout << "Value[" << index << "]: " << value << endl;
          index++;
      }

      // Reset stream to process it for value assignment
      lineStream.clear(); // Clear any error flags
      lineStream.seekg(0, ios::beg); // Reset position
    */

    // Set stem migration parameters m1-m2
    for (int i = 0; i < 2; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Stem::CaAlgMigration[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Stem::CaAlgMigration[" << i << "] = " << Stem::CaAlgMigration[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to stem migration parameter #" << i << endl;
        }
    }

    // Set stem proliferation parameters m3-m6
    for (int i = 0; i < 4; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Stem::proliferation[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Stem::proliferation[" << i << "] = " << Stem::proliferation[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to stem proliferation parameter #" << i << endl;
        }
    }

    // Set stem cytokine synthesis parameters m7-m9
    for (int i = 0; i < 3; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Stem::cytokineSynthesis[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Stem::cytokineSynthesis[" << i << "] = " << Stem::cytokineSynthesis[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to stem cytokine synthesis parameter #" << i << endl;
        }
    }

    // Set stem collagen synthesis parameters m10
    if (getline(lineStream, value, '\n')) {
        float value_as_float = atof(value.c_str());
        Stem::CollagenSynth[0] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
        cout << "Stem::CollagenSynth" << " = " << Stem::CollagenSynth[0] << endl;
#endif
    }
    else {
        cerr << "Error in assigning value to stem collagen synthesis parameter" << endl;
    }

    // Set stem aggrecan synthesis parameters m11
    if (getline(lineStream, value, '\n')) {
        float value_as_float = atof(value.c_str());
        Stem::AggrecanSynth[0] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
        cout << "Stem::AggrecanSynth" << " = " << Stem::AggrecanSynth[0] << endl;
#endif
    }
    else {
        cerr << "Error in assigning value to stem aggrecan synthesis parameter" << endl;
    }

    // Set stem differentiation parameters m12-m16
    for (int i = 0; i < 5; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Stem::differentiation[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Stem::differentiation[" << i << "] = " << Stem::differentiation[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to stem differentiation parameter #" << i << endl;
        }
    }

    // Set pre-NP migration parameters p1-p2
    for (int i = 0; i < 2; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Progen::CaAlgMigration[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Progen::CaAlgMigration[" << i << "] = " << Progen::CaAlgMigration[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to pre-np migration parameter #" << i << endl;
        }
    }

    // Set pre-NP proliferation parameter p3
    if (getline(lineStream, value, '\n')) {
        float value_as_float = atof(value.c_str());
        Progen::proliferation[0] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
        cout << "Progen::proliferation" << " = " << Progen::proliferation[0] << endl;
#endif
    }
    else {
        cerr << "Error in assigning value to pre-NP proliferation parameter" << endl;
    }

    // Set pre-NP cytokine synthesis parameters p4-p6
    for (int i = 0; i < 3; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Progen::cytokineSynthesis[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Progen::cytokineSynthesis[" << i << "] = " << Progen::cytokineSynthesis[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to pre-np cytokine synthesis parameter #" << i << endl;
        }
    }

    // Set pre-NP aggrecan synthesis parameter p7
    if (getline(lineStream, value, '\n')) {
        float value_as_float = atof(value.c_str());
        Progen::AggrecanSynth[0] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
        cout << "Progen::AggrecanSynth" << " = " << Progen::AggrecanSynth[0] << endl;
#endif
    }
    else {
        cerr << "Error in assigning value to pre-NP aggrecan synthesis parameter" << endl;
    }

    // Set pre-NP differentiation parameters p8-p10
    for (int i = 0; i < 3; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Progen::differentiation[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Progen::differentiation[" << i << "] = " << Progen::differentiation[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to pre-np differentiation parameter #" << i << endl;
        }
    }

    // Set NP migration parameters k0-k1
    for (int i = 0; i < 2; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            NP::CaAlgMigration[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "NP::CaAlgMigration[" << i << "] = " << NP::CaAlgMigration[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to NP migration parameter #" << i << endl;
        }
    }

    // Set NP proliferation parameters k2-k7
    for (int i = 0; i < 6; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Cell::proliferation[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Cell::proliferation[" << i << "] = " << Cell::proliferation[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to NP proliferation parameter #" << i << endl;
        }
    }

    // Set NP/general cytokine synthesis parameters k8-k17
    for (int i = 0; i < 10; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            Cell::cytokineSynthesis[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "Cell::cytokineSynthesis[" << i << "] = " << Cell::cytokineSynthesis[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to NP cytokine synthesis parameter #" << i << endl;
        }
    }

    // Set NP collagen synthesis parameters k18-k20
    for (int i = 0; i < 3; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            NP::CollagenSynth[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "NP::CollagenSynth[" << i << "] = " << NP::CollagenSynth[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to NP collagen synthesis parameter #" << i << endl;
        }
    }

    // Set NP aggrecan synthesis parameters k21-k23
    for (int i = 0; i < 3; ++i) {
        if (getline(lineStream, value, '\n')) {
            float value_as_float = atof(value.c_str());
            NP::AggrecanSynth[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
            cout << "NP::AggrecanSynth[" << i << "] = " << NP::AggrecanSynth[i] << endl;
#endif
        }
        else {
            cerr << "Error in assigning value to NP aggrecan synthesis parameter #" << i << endl;
        }
    }

  /* // Set WHWorld cytokine decay parameters
    for (int i = 0; i < 8; ++i) {
     if((i!=3) && (i!=4) && (i!=7)){
       if (getline(lineStream, value, '\t')) {
          float value_as_float = atof(value.c_str());
          
          // Convert half-life to decay rate
          if (value_as_float > 0) {
          float decay_rate = pow(0.5, 30/value_as_float);
          WHWorld::cytokineDecay[i] = decay_rate;
          WHWorld::halfLifes_static[i] = value_as_float;
          #ifdef PRINT_PARAMETER_VALUES
            cout << "WHWorld::halfLifes_static[" << i << "] = " << value_as_float << ", WHWorld::cytokineDecay["<< i << "] = " << WHWorld::cytokineDecay[i] << endl;
          #endif
         } else {
           cerr << "Error in assigning value to WHWorld half life parameter and cytokine decay parameter #" << i << " (Half-life <= 0)" << endl;
         }
       } else {
         cerr << "Error in assigning value to WHWorld half life parameter and cytokine decay parameter #" << i << endl;
       }
     }
   } 
  */

  // Set WHWorld elastic modulus parameters c1-c7
  for (int i = 0; i < 7; ++i) {
    if (getline(lineStream, value, '\n')) {
      float value_as_float = atof(value.c_str());
     WHWorld::ElasticMod[i] = value_as_float;
    #ifdef PRINT_PARAMETER_VALUES
      cout << "WHWorld::ElasticMod[" << i << "] = " << WHWorld::ElasticMod[i] << endl;
    #endif
    } else {
      cerr << "Error in assigning value to elastic modulus parameter #"<< i << endl;
    }
  }

  // Set WHWorld pore size parameters c8-c9
  for (int i = 0; i < 2; ++i) {
      if (getline(lineStream, value, '\n')) {
          float value_as_float = atof(value.c_str());
          WHWorld::PoreSize[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
          cout << "WHWorld::PoreSize[" << i << "] = " << WHWorld::PoreSize[i] << endl;
#endif
      }
      else {
          cerr << "Error in assigning value to pore size parameter #" << i << endl;
      }
  }

  // Set WHWorld mass loss parameters c10-c13
  for (int i = 0; i < 4; ++i) {
      if (getline(lineStream, value, '\n')) {
          float value_as_float = atof(value.c_str());
          WHWorld::MassLoss[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
          cout << "WHWorld::MassLoss[" << i << "] = " << WHWorld::MassLoss[i] << endl;
#endif
      }
      else {
          cerr << "Error in assigning value to mass loss parameter #" << i << endl;
      }
  }

  // Set WHWorld swelling ratio parameters c14-c18
  for (int i = 0; i < 5; ++i) {
      if (getline(lineStream, value, '\n')) {
          float value_as_float = atof(value.c_str());
          WHWorld::SwellRatio[i] = value_as_float;
#ifdef PRINT_PARAMETER_VALUES
          cout << "WHWorld::SwellRatio[" << i << "] = " << WHWorld::SwellRatio[i] << endl;
#endif
      }
      else {
          cerr << "Error in assigning value to swelling ratio parameter #" << i << endl;
      }
  }

/*
  // Set WHWorld XL density parameters
  for (int i = 0; i < 2; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      WHWorld::XLDensity[i] = value_as_float;
      #ifdef PRINT_PARAMETER_VALUES
        cout << "WHWorld::XLDensity[" << i << "] = " << WHWorld::XLDensity[i] << endl;
      #endif
    } else {
      cerr << "Error in assigning value to XL density parameter #"<< i << endl;
    }
  }
*/

  /* // Set Agent migration parameters
    for (int i = 0; i < 2; ++i) {
      if (getline(lineStream, value, '\t')) {
        float value_as_float = atof(value.c_str());
        Agent::CaAlgMigration[i] = value_as_float;
        #ifdef PRINT_PARAMETER_VALUES
          cout << "Agent::CaAlgMigration[" << i << "] = " << Agent::CaAlgMigration[i] << endl;
        #endif
      } else {
        cerr << "Error in assigning value to agent Ca-Alg migration parameter #"<< i << endl;
      }
    } 
  */

/*
  // Set Agent CaAlg proliferation parameters 
  for (int i = 0; i < 5; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Agent::CaAlgProlif[i] = value_as_float;
      #ifdef PRINT_PARAMETER_VALUES
        cout << "Agent::CaAlgProlif[" << i << "] = " << Agent::CaAlgProlif[i] << endl;
      #endif
    } else {
      cerr << "Error in assigning value to agent CaAlg proliferation parameter #"<< i << endl;
    }
  }
*/

/*
  // Set Agent CaAlg viability parameters 
  for (int i = 0; i < 3; ++i) {
    if (getline(lineStream, value, '\t')) {
      float value_as_float = atof(value.c_str());
      Agent::CaAlgViability[i] = value_as_float;
      #ifdef PRINT_PARAMETER_VALUES
        cout << "Agent::CaAlgViability[" << i << "] = " << Agent::CaAlgViability[i] << endl;
      #endif
    } else {
      cerr << "Error in assigning value to agent CaAlg viability parameter #"<< i << endl;
    }
  }
*/

  /* // Set Agent viability rate parameters 
    for (int i = 0; i < 1; ++i) {
      if (getline(lineStream, value, '\t')) {
       float value_as_float = atof(value.c_str());
       Agent::viabilityRate = value_as_float;
       #ifdef PRINT_PARAMETER_VALUES
         cout << "Agent::viabilityRate = " << Agent::viabilityRate << endl;
       #endif
     } else {
       cerr << "Error in assigning value to agent viability rate parameter #"<< i << endl;
     }
   } 
  */

/* // Set Agent proliferation rate parameters 
   for (int i = 0; i < 1; ++i) {
    if (getline(lineStream, value, '\t')) {
       float value_as_float = atof(value.c_str());
       Agent::proliferationRate = value_as_float;
        #ifdef PRINT_PARAMETER_VALUES
              cout << "Agent::proliferationRate = " << Agent::proliferationRate << endl;
        #endif
     } else {
       cerr << "Error in assigning value to agent proliferation rate parameter #"<< i << endl;
    }
   } 
*/

input_file.close();
#endif  // ifdef CALIBRATION

return;
}

void outputTotalChem(WHWorld* myWorld, string filename) {
#ifdef CALIBRATION
  cout << "Outputting total chem to: " << filename << endl;
	//ofstream output_file(filename, ios::app);
  ofstream output_file;
  output_file.open(filename.c_str(), ios::app);

  //Get collagen & total protein amount:
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

  tp += myWorld->WHWorldChem.totalTNF;
  tp += myWorld->WHWorldChem.totalTGF;
  tp += myWorld->WHWorldChem.totalIL1beta;

  //Output:
  output_file << fixed << myWorld->cells.actualSize() << "\t";
  output_file << col << "\t";
  output_file << agg << "\t";
  output_file << tp << "\t";

  output_file << fixed << myWorld->WHWorldChem.totalTNF << "\t";
  output_file << myWorld->WHWorldChem.totalTGF << "\t";
  output_file << myWorld->WHWorldChem.totalIL1beta << "\t";
  
  output_file << endl;
  
#endif  // ifdef CALIBRATION
	return;
}
}  // namespace util

#endif /* PARAMETERS_H_ */