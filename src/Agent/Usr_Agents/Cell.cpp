/* 
 * File: Cell.cpp
 *
 * File Contents: Contains the Cell class.
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *               Meghana Munipalle
 */

#include "Cell.h"
#include "../../enums.h"
#include <iostream>
#include <algorithm>      
#include <cmath>  

using namespace std;
int Cell::numOfCells = 0;
float Cell::proliferation[6] = { 24, 10, 1, 0, 25, 2 };
float Cell::cytokineSynthesis[10] = {10, 0.05, 10, 5, 2.4, 4, 2, 5, 1, 3.2}; //calibration variables
float Cell::activation[5] = {500, 50, 0, 25, 2.5}; // these are in sample.txt, but unutilized in new version
float Cell::ECMsynthesis[12] = { 1, 1, 1, 50, 25, 2, 10, 5, 1, 1, 25, 2 }; // these are in sample.txt, but unutilized in new version

int Stem::numOfStem = 0;
float Stem::migrationSpeed = 1; // patch/tick
float Stem::apoptosisChance = 0.05;
float Stem::collagenSynthRate = 1; // placeholder values which will be recalculated
float Stem::aggrecanSynthRate = 0.5; // placeholder values which will be recalculated

float Stem::CaAlgMigration[2] = { 0.11, 0.35 };
float Stem::cytokineSynthesis[3] = { 50, 0, -0.807 };
//float Stem::ECMsynthesis[4] = {};
float Stem::CollagenSynth[1] = { 10 };
float Stem::AggrecanSynth[1] = { 100000 };
float Stem::proliferation[4] = {10, 24, 0.8, 0.001};
float Stem::differentiation[5] = { 0.7, 0.3, 0.5, 0.001, 48 };

int Progen::numOfProgen = 0; 
float Progen::migrationSpeed = 1;    // patch/tick
float Progen::apoptosisChance = 0.1;
float Progen::aggrecanSynthRate = 1;

float Progen::CaAlgMigration[2] = { 0.11, 0.83 };
float Progen::cytokineSynthesis[3] = { 50, 0, -0.807 };
float Progen::AggrecanSynth[1] = { 1 };
float Progen::proliferation[1] = {24};
float Progen::differentiation[3] = {0.7, 0.3, 48};

int NP::numOfNP = 0;
float NP::migrationSpeed = 1;    // patch/tick
float NP::collagenSynthRate = 1;
float NP::aggrecanSynthRate = 1.5;

float NP::CaAlgMigration[2] = { 0.11, 1.30 };
float NP::CollagenSynth[3] = { 10, 6.45, 3.6 };
float NP::AggrecanSynth[3] = { 20, 38, 16.6 };


Cell::Cell() {
	cout << "default cell alloc" << endl;

	// added for debugging: print out what cell type
	if (typeid(*this) == typeid(Stem)) {
		cout << "cell type stem";
	}
	else if (typeid(*this) == typeid(Progen)) {
		cout << "cell type progen";
	}
	else if (typeid(*this) == typeid(NP)) {
		cout << "cell type NP";
	}
}

Cell::Cell(Patch* patchPtr) {
	this->ix[write_t] = patchPtr->indice[0];
	this->iy[write_t] = patchPtr->indice[1];
	this->iz[write_t] = patchPtr->indice[2];
	this->index[write_t] = patchPtr->index;
	this->alive[write_t] = true;

	this->activate[write_t] = false;
	//this->color[write_t] = ccell;
	this->size[write_t] = 2;
	//this->type[write_t] = cell;
	this->ix[read_t] = patchPtr->indice[0];
	this->iy[read_t] = patchPtr->indice[1];
	this->iz[read_t] = patchPtr->indice[2];
	this->index[read_t] = patchPtr->index;

	/* In OMP version, we wait to add cells at the end of the tick, whereas in serial version, cells are always added right away.
     * Thus, cells, when added in OMP, should be alive right away. */
	#ifdef _OMP
		this->alive[read_t] = true;
		this->life[read_t] = this->life[write_t];
	#else
		this->alive[read_t] = false;
		this->life[read_t] = 0;
	#endif
	
	this->activate[read_t] = false;
	//this->color[read_t] = ccell;
	this->size[read_t] = 2;
	//this->type[read_t] = cell;
	/* Added by MM to check types of cell stages and add to respective counters: */
	if (typeid(*this) == typeid(Stem)) {
		Stem::numOfStem++;
	} else if (typeid(*this) == typeid(Progen)) {
		Progen::numOfProgen++;
	} else if (typeid(*this) == typeid(NP)) {
		NP::numOfNP++;
	}
	Cell::numOfCells++;
}

Stem::Stem(Patch* patchPtr) {
	//Cell::Cell(patchPtr);

	this->color[write_t] = cstem;
	this->type[write_t] = stem;

	this->color[read_t] = cstem;
	this->type[read_t] = stem;
}

Progen::Progen(Patch* patchPtr) {
	//Cell::Cell(patchPtr);

	this->color[write_t] = cprogen;
	this->type[write_t] = progen;

	this->color[read_t] = cprogen;
	this->type[read_t] = progen;
}

NP::NP(Patch* patchPtr) {
	//Cell::Cell(patchPtr);

	this->color[write_t] = cnp;
	this->type[write_t] = np;

	this->color[read_t] = cnp;
	this->type[read_t] = np;

	#ifndef MODEL_SCAFFOLD
		// Unactivated chondrocytes live for 5 to 11 days. 0 corresponds to hours:
		if (Agent::agentWorldPtr->clock == 0) this->life[write_t] = WHWorld::reportTick(0, rand() % 12);
		else this->life[write_t] = WHWorld::reportTick(0, 5 + rand() % 7);
	#else
		this->life[write_t] = WHWorld::reportTick(0, 5 + rand() % 7);
	#endif
}

Cell::Cell(int x, int y, int z) {
	this->ix[write_t] = x;
	this->iy[write_t] = y;
	this->iz[write_t] = z;
	this->index[write_t] = x + y*nx + z*nx*ny;
	this->alive[write_t] = true;

	#ifndef MODEL_SCAFFOLD
		// Unactivated chondrocytes live for 5 to 11 days. 0 corresponds to hours.
		if (Agent::agentWorldPtr->clock == 0) this->life[write_t] = WHWorld::reportTick(0, rand()%12);
		else this->life[write_t] = WHWorld::reportTick(0, 5 + rand()%7);
	#else
		this->life[write_t] = WHWorld::reportTick(0, 5 + rand()%7);
	#endif

	this->activate[write_t] = false;
	this->color[write_t] = ccell;
	this->size[write_t] = 2;
	this->type[write_t] = cell;

  /* In OMP version, we wait to add cells at the end of the tick, whereas in serial version, cells are always added right away.
   * Thus, cells, when added in OMP, should be alive right away. */
	#ifdef _OMP
		this->ix[read_t] = x;
		this->iy[read_t] = y;
		this->iz[read_t] = z;
		this->index[read_t] = x + y*nx + z*nx*ny;
		this->alive[read_t] = true;
		this->life[read_t] = this->life[write_t];
		this->activate[read_t] = false;
		this->color[read_t] = ccell;
		this->size[read_t] = 2;
		this->type[read_t] = cell;
	#else
		this->ix[read_t] = x;
		this->iy[read_t] = y;
		this->iz[read_t] = z;
		this->index[read_t] = x + y*nx + z*nx*ny;
		this->alive[read_t] = false;
		this->life[read_t] = 0;
		this->activate[read_t] = false;
		this->color[read_t] = ccell;
		this->size[read_t] = 2;
		this->type[read_t] = cell;
	#endif
	/* Added by MM to check types of cell stages and add to respective counters: */
	if (typeid(*this) == typeid(Stem)) {
		Stem::numOfStem++;
	}
	else if (typeid(*this) == typeid(Progen)) {
		Progen::numOfProgen++;
	}
	else if (typeid(*this) == typeid(NP)) {
		NP::numOfNP++;
	}
	Cell::numOfCells++;  
}

Stem::Stem(int x, int y, int z) {}

Progen::Progen(int x, int y, int z) {}

NP::NP(int x, int y, int z) {}

Cell::~Cell() {}

void Cell::cellFunction() {
	if (this->alive[read_t] == false) return;
	if (this->alive[read_t] == true) {
		//Cell temp;
		Cell* temp = this;
		if (typeid(*this) == typeid(Stem)) {
			Stem* tmp = dynamic_cast<Stem*>(temp); // dynamic cast to next derived class (cell type)
			//Stem tmp = dynamic_cast<Stem*>(temp); // dynamic cast to next derived class (cell type)
			if (tmp == nullptr) {
				cout << "Casting Failed" << endl;
			}
			tmp->Stem::stem_cellFunction();
		} else if (typeid(*this) == typeid(Progen)) {
			Progen* tmp = dynamic_cast<Progen*>(temp); // dynamic cast to next derived class (cell type)
			if (tmp == nullptr) {
				cout << "Casting Failed" << endl;
			}
			tmp->Progen::progen_cellFunction();
		} else if (typeid(*this) == typeid(NP)) {
			NP* tmp = dynamic_cast<NP*>(temp); // dynamic cast to next derived class (cell type)
			if (tmp == nullptr) {
				cout << "Casting Failed" << endl;
			}
			tmp->NP::NP_cellFunction();
		}
	}
	//if (this->activate[read_t] == false) this->chond_cellFunction();
	//else this->achond_cellFunction();
}

void NP::NP_cellFunction() {

	int in = this->index[read_t];
	double hours = Agent::agentWorldPtr->reportHour();
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

	/* Unactivated chondrocytes only move along their preferred gradient and perform biological functions if there is damage. */
	if (totaldamage == 0) {

#ifdef MODEL_SCAFFOLD
		// If cell is moving on Ca-Alg substrate chondrocyte move up to "migrationSpeed" patches per tick determined by substrate composition
		if (NP::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
			for (int dx = 0; dx < NP::migrationSpeed; dx++) this->wiggle();

		}
		else if (rollDice(0.25)) { // If cell is not actively migrating, consider chance of moving to next patch            
			this->wiggle();
		}
#else
		this->wiggle();
#endif

	}
	else {

		/* -------------------------------------------------------------------------- */
		/*                                PROLIFERATION                               */
		/* -------------------------------------------------------------------------- */
#ifdef MODEL_SCAFFOLD
	// Cells in CaAlg hydrogel proliferate at a proliferation rate given time and hydrogel composition (% Alg) 
		in = this->index[read_t];
		if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) {

#ifdef CALIBRATION
			if ((fmod((float)Agent::agentWorldPtr->clock, 2) == 0) && rollDice(Agent::proliferationRate)) {   //check every hour 
#else 
			if ((fmod((float)Agent::agentWorldPtr->clock, 2) == 0) && rollDice(Agent::proliferationRate)) { //check every hour
#endif

				float meanTNF = this->meanNeighborChem(TNF);
				float meanTGF = this->meanNeighborChem(TGF);
				float meanIL1 = this->meanNeighborChem(IL1beta);
				int countfHA = 0;
				int TGFrelated = 0;

#ifndef CALIBRATION
				if (meanTGF <= Cell::proliferation[1]) {
#else 
				if (meanTGF <= 10) {
#endif
					TGFrelated = 1;  // Low TGF (0.1-1 ng) stimulate chond proliferation and attraction. 
				}
				else {
					TGFrelated = -1; // High TGF (1-10 ng) inhibits proliferation.
				}

#ifndef CALIBRATION
				float cellProlif = Cell::proliferation[2] * (log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF)) + Cell::proliferation[3];
				if (rollDice(Cell::proliferation[4] + cellProlif / Cell::proliferation[5])) {
#else
				float cellProlif = 1 * log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF) + 0;
				if (rollDice(25 + cellProlif / 2)) {
#endif	
					this->Agent::hatchnewcell(1, np);
					this->die();
					return;
				}
			}
			// Unactivated chondrocytes proliferate every 24 hours.
		}
		else {
#endif

#ifndef CALIBRATION
			if (fmod((float)hours, Cell::proliferation[0]) == 0) {
#else  
			if (fmod(hours, 24) == 0) {
#endif  
				float meanTNF = this->meanNeighborChem(TNF);
				float meanTGF = this->meanNeighborChem(TGF);
				float meanIL1 = this->meanNeighborChem(IL1beta);
				int countfHA = this->countNeighborECM(fha);
				int TGFrelated = 0;

#ifndef CALIBRATION
				if (meanTGF <= Cell::proliferation[1]) {
#else
				if (meanTGF <= 10) {
#endif
					TGFrelated = 1;  // Low TGF (0.1-1 ng) stimulate chond proliferation and attraction.
				}
				else {
					TGFrelated = -1;  // High TGF (1-10 ng) inhibits proliferation.
				}

#ifndef CALIBRATION
				float cellProlif = Cell::proliferation[2] * (log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF)) + Cell::proliferation[3];
				if (rollDice(Cell::proliferation[4] + cellProlif / Cell::proliferation[5])) {
#else  
				float cellProlif = log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF);
				if (rollDice(25 + cellProlif / 2)) {
#endif 
					this->Agent::hatchnewcell(1, np);
					this->die();
					return;
				}
			}
#ifdef MODEL_SCAFFOLD
		}
#endif

		/* -------------------------------------------------------------------------- */
		/*                                  MOVEMENT                                  */
		/* -------------------------------------------------------------------------- */
		this->cellSniff();

		/* -------------------------------------------------------------------------- */
		/*                                 ACTIVATION                                 */
		/* -------------------------------------------------------------------------- */
		//// An unactivated chondrocyte can be activated if it is in the damage zone:
		//if (Agent::agentPatchPtr[in].inDamzone == true) {
		//	// Low TGF promote and high TGF inhibit chances of chondrocyte activation:
		//	int patchTGF = agentWorldPtr->WHWorldChem.pTGF[in];
		//	
		//#ifndef CALIBRATION
		//	if ((patchTGF > 10 && rollDice(Cell::activation[1])) || (patchTGF > Cell::activation[2]) || (rollDice(Cell::activation[3]))) { 
		//#else  
		//	if ((patchTGF > 10 && rollDice(50.0)) || (patchTGF > 0) || (rollDice(25))) { 
		//#endif
		//		this->chondActivation();
		//	}
		//	}
		//}

		/* -------------------------------------------------------------------------- */
		/*                      ECM PROTEIN & CHEMICAL SYNTHESIS                      */
		/* -------------------------------------------------------------------------- */
		// Calculates chemical gradients and patch chemical concentrations:
		//	int in = this->index[read_t];
		float meanTNF = this->meanNeighborChem(pTNF);
		float meanTGF = this->meanNeighborChem(pTGF);
		float meanIL1 = this->meanNeighborChem(pIL1beta);
		int countnHA = 0;
		int countfHA = 0;
		float patchTNF = this->agentWorldPtr->WHWorldChem.pTNF[in];
		float patchTGF = this->agentWorldPtr->WHWorldChem.pTGF[in];
		float patchIL1beta = (this->agentWorldPtr->WHWorldChem.pIL1beta[in]);

		// Makes collagen and aggrecan every 12 hours.
#ifndef CALIBRATION
		if (fmod(((Agent::agentWorldPtr)->reportHour()), 1) == 0) {
#else 
		if (fmod(((Agent::agentWorldPtr)->reportHour()), 1) == 0) { //12
#endif

#ifdef MODEL_SCAFFOLD
			// Active cell adhered to Ca-Alg synthesis ECM according the substrate mechanical properties:
			int in = this->index[read_t];
			if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
				for (int i = 0; i < NP::collagenSynthRate; i++) {
					this->makeOCollagen(meanTGF, meanIL1);
				}

				for (int i = 0; i < NP::aggrecanSynthRate; i++) {
					this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
				}

			}
			else {
				this->makeOCollagen(meanTGF, meanIL1);
				this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
			}

#else
			this->makeOCollagen(meanTGF, meanIL1);
			this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
#endif
		}

		// Change in chemicals due to cells:
#ifndef CALIBRATION
		(this->agentWorldPtr->WHWorldChem.dTGF[in]) += Cell::cytokineSynthesis[0] + Cell::cytokineSynthesis[1] * (patchTGF) + Cell::cytokineSynthesis[2]*(patchIL1beta) + Cell::cytokineSynthesis[3]*(patchTNF);			//(this->agentWorldPtr->WHWorldChem.dTGF[in]) +=  Chondrocyte::cytokineSynthesis[0] + Chondrocyte::cytokineSynthesis[1]*(1 + Chondrocyte::cytokineSynthesis[2]*patchTNF);
		(this->agentWorldPtr->WHWorldChem.dTNF[in]) += Cell::cytokineSynthesis[4] + (Cell::cytokineSynthesis[5] * ((patchIL1beta) / (1 + Cell::cytokineSynthesis[6] * patchTGF));	//(this->agentWorldPtr->WHWorldChem.dTNF[in]) += Chondrocyte::cytokineSynthesis[3] + Chondrocyte::cytokineSynthesis[4]/(1 + patchTGF*Chondrocyte::cytokineSynthesis[5]);
		(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += Cell::cytokineSynthesis[7] + (Cell::cytokineSynthesis[8] * ((patchTNF) / (1 + Cell::cytokineSynthesis[9] * patchTGF)); //(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += Chondrocyte::cytokineSynthesis[6] + (Chondrocyte::cytokineSynthesis[7]*patchTNF)/(Chondrocyte::cytokineSynthesis[8] + Chondrocyte::cytokineSynthesis[9]*patchTGF);
#else
		(this->agentWorldPtr->WHWorldChem.dTGF[in]) += 10 + 0.05 * (patchTGF + 10 * patchTNF); 				//9.98 + 2.58*patchTGF + 5.11*patchTNF;				//2.11 + 3.7*patchTGF;
		(this->agentWorldPtr->WHWorldChem.dTNF[in]) += 5 + (2.4 * patchIL1beta) / (1 + 4 * patchTGF);				//5.16 + (2.42*patchIL1beta)/(1 + 4.22*patchTGF);	//2.4*patchIL1beta + 4.8/(1 + 1.27*patchTGF);		
		(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += 2 + (5 * patchTNF) / (1 + 3.2 * patchTGF);		//2.11 + (5.43*patchTNF)/(1 + 3.26*patchTGF);		//4;
#endif

		/* -------------------------------------------------------------------------- */
		/*                                    DEATH                                   */
		/* -------------------------------------------------------------------------- */

		// Cells in Ca-Alg hydrogel have at a viability/death rate determined by time:
#ifdef MODEL_SCAFFOLD
#ifdef CALIBRATION
		if (fmod((float)Agent::agentWorldPtr->clock, Agent::CaAlgViability[2]) == 0 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
#else
		if (fmod((float)Agent::agentWorldPtr->clock, 6) && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
#endif
			if (rollDice(100 - Agent::viabilityRate)) {
				this->die();
				return;
			}
		}
#endif

		// Unactivated chondrocytes can die naturally:
		this->life[write_t] = this->life[read_t] - 1;
		if (this->life[read_t] <= 0) this->die();
	}
}

void Cell::cellSniff() {
	int in = this->index[read_t];
	if ((Agent::agentPatchPtr[in]).inDamzone == true) {
		if (rollDice(80) && this->moveToHighestChem(pcellgrad) == true){
			#ifdef MODEL_SCAFFOLD
			if (typeid(*this) == typeid(Stem)) {
				// If cell is moving on Ca-Alg substrate:
				if (Stem::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
					// Move up to "migrationSpeed" patches per tick:
					for (int dx = 0; dx < Stem::migrationSpeed; dx++) this->wiggle();
				}
			}
			else if (typeid(*this) == typeid(Progen)) {
				if (Progen::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
					// Move up to "migrationSpeed" patches per tick:
					for (int dx = 0; dx < Progen::migrationSpeed; dx++) this->wiggle();
				}
			}
			else if (typeid(*this) == typeid(NP)) {
				if (NP::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
					// Move up to "migrationSpeed" patches per tick:
					for (int dx = 0; dx < NP::migrationSpeed; dx++) this->wiggle();
				}
			}
				// If cell is not actively migrating, consider chance of moving to next patch:
		} else if(rollDice(0.25)){		
			this->wiggle(); 
		}
		#else
			this->wiggle();
		#endif

	} else {
		#ifdef MODEL_SCAFFOLD
			// If cell is moving on Ca-Alg substrate:
		if (typeid(*this) == typeid(Stem)) {
			// If cell is moving on Ca-Alg substrate:
			if (Stem::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
				// Move up to "migrationSpeed" patches per tick:
				for (int dx = 0; dx < Stem::migrationSpeed; dx++) this->wiggle();
			}
		}
		else if (typeid(*this) == typeid(Progen)) {
			if (Progen::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
				// Move up to "migrationSpeed" patches per tick:
				for (int dx = 0; dx < Progen::migrationSpeed; dx++) this->wiggle();
			}
		}
		else if (typeid(*this) == typeid(NP)) {
			if (NP::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
				// Move up to "migrationSpeed" patches per tick:
				for (int dx = 0; dx < NP::migrationSpeed; dx++) this->wiggle();
			}
		}
		else if (rollDice(0.25)) this->wiggle(); 							   // If cell is not actively migrating, consider chance of moving to next patch
		#else
			this->wiggle();
		#endif
	}

	// TGF can excite NP cell and overcome gradient: //NOTE MM: double check this for MSCs
	if ((typeid(*this) == typeid(NP) && this->meanNeighborChem(pTGF) > 0)) {
		#ifdef MODEL_SCAFFOLD
			// If cell is moving on Ca-Alg substrate:
			if (NP::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
				for (int dx = 0; dx < NP::migrationSpeed; dx++) this->wiggle(); // Move up to "migrationSpeed" patches per tick:
			}
			else if (rollDice(0.25)) { this->wiggle(); } 						   // If cell is not actively migrating, consider chance of moving to next patch: 		
		#else
			this->wiggle();
		#endif 
	}
}

void Cell::die() {
	int in = this->index[read_t];
	Agent::agentPatchPtr[in].clearOccupied();
	Agent::agentPatchPtr[in].occupiedby[write_t] = nothing;
	this->alive[write_t] = false;
	this->life[write_t] = 0;
}

//void Chondrocyte::chondActivation() {
//	int in = this->index[read_t];
//	int target = this->index[write_t]; // This assumes that after this function is called, no more move() would be called in the same tick and the cell will not naturally die
//	if (this->activate[read_t] == false && this->life[read_t] > 1) {
//		this->type[write_t] = achondrocyte;
//		this->activate[write_t] = true;
//		this->color[write_t] = cachondrocyte;
//		Agent::agentPatchPtr[target].setOccupied();
//		Agent::agentPatchPtr[target].occupiedby[write_t] = achondrocyte;
//	}
//}

//void Cell::cellDeactivation() {
//	int in = this->index[read_t];
//	int target = this->index[write_t]; // This assumes that after this function is called, no more move() would be called in the same tick 
//	if (this->activate[read_t] == true&& this->life[read_t] > 1) {
//		this->type[write_t] = stem;
//		this->activate[write_t] = false;
//		this->color[write_t] = cstem;
//		Agent::agentPatchPtr[in].setOccupied();
//		Agent::agentPatchPtr[in].occupiedby[write_t] = stem;
//	}
//}

void Cell::copyAndInitialize(Agent* original, int dx, int dy, int dz) {
	int in = this->index[read_t];

	// Initializes location of new Cell relative to original agent:
	this->ix[write_t] = original->getX() + dx;
	this->iy[write_t] = original->getY() + dy;
	this->iz[write_t] = original->getZ() + dz;
	this->index[write_t] = this->ix[write_t] + this->iy[write_t]*Agent::nx + this->iz[write_t]*Agent::nx*Agent::ny;
  	
	// Initializes new Cell:
	this->alive[read_t] = true;
	this->life[read_t] = WHWorld::reportTick(0, 5 + rand()%7);  // Unactivated chondrocytes live for 5 to 11 days. 0 corresponds to hours.
	this->activate[read_t] = false;
	this->color[read_t]= ccell;
	this->size[read_t] = 2;
	this->type[read_t] = cell;
	this->alive[write_t] = true;
	this->life[write_t] = this->life[read_t];
	this->activate[write_t] = false;
	this->color[write_t]= ccell;
	this->size[write_t] = 2;
	this->type[write_t] = cell;

	Cell::numOfCells++;

  	// Assigns new Cell to this patch if it is unoccupied:
	if (Agent::agentPatchPtr[in].isOccupied() == false) {
		Agent::agentPatchPtr[in].setOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = this->type[read_t];
	} else {
		cout << "error in hatching and initialization!!!" << dx << " " << dy << endl;
	}
}

void Stem::stem_cellFunction() {
	int in = this->index[read_t];
	double hours = Agent::agentWorldPtr->reportHour();
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

	/* -------------------------------------------------------------------------- */
	/*                                PROLIFERATION                               */
	/* -------------------------------------------------------------------------- */
  	
	// Stem cells in vitro (MODEL_SCAFFOLD) proliferate:
	#ifdef MODEL_SCAFFOLD
		in = this->index[read_t];
		cout << "agent patch type is " << Agent::agentPatchPtr[in].type[read_t] << endl; //added for debug
		if (Agent::agentPatchPtr[in].type[read_t] == CaAlg){
		#ifdef CALIBRATION
			if (fmod((float)hours, Stem::proliferation[1]) == 0) {
		#else 
			if (fmod(hours, 24) == 0) {
		#endif 
				float meanTNF = this->meanNeighborChem(TNF);
				float meanTGF = this->meanNeighborChem(TGF);
				float meanIL1 = this->meanNeighborChem(IL1beta);
				//int countfHA = this->countNeighborECM(fha);
				int TGFrelated = 0;
				#ifndef CALIBRATION
					if (meanTGF <= Stem::proliferation[0]) {
				#else  
					if (meanTGF <= 10) {
				#endif 

						TGFrelated = 1;  // Low TGF (0.1-1ng) stimulate proliferation and attraction. 
					} else {
						TGFrelated = -1; // High TGF (1-10ng) inhibits proliferation. 
					}

					#ifndef CALIBRATION
						float stemProlif = log10(1 - Stem::proliferation[2]*meanTNF - Stem::proliferation[3]*meanIL1 + TGFrelated*meanTGF);
						if (rollDice(stemProlif) {  
					#else  
						float stemProlif = log10(1 + meanTNF + meanIL1 + TGFrelated*meanTGF); 
						if (rollDice(stemProlif)) { 										
					#endif  
							this->Agent::hatchnewcell(1, stem);
							//this->die();
							return;
						}
						}
					} else {
	#endif

			#ifndef CALIBRATION
				if (fmod((float)hours, Stem::proliferation[1])) == 0 {
			#else  
				if (fmod(hours, 24) == 0) {
			#endif  
					float meanTNF = this->meanNeighborChem(TNF);
					float meanTGF = this->meanNeighborChem(TGF);
					float meanIL1 = this->meanNeighborChem(IL1beta);
					//int countfHA = this->countNeighborECM(fha);
					int TGFrelated = 0;

  				#ifndef CALIBRATION
					if (meanTGF <= Stem::proliferation[0]) {
  				#else  
					if (meanTGF <= 10) {
  				#endif  
						TGFrelated = 1; // Low TGF (0.1-1nm) stimulate chond proliferation and attraction. 
					} else {
						TGFrelated = -1;  // High TGF (1-10nm) inhibits proliferation. 
					}

				#ifndef CALIBRATION
					float stemProlif = log10(1 - Stem::proliferation[2] * meanTNF - Stem::proliferation[3] * meanIL1 + TGFrelated * meanTGF);
					if (rollDice(stemProlif) {
				#else
					float stemProlif = log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF);
						if (rollDice(stemProlif)) {
				#endif
						this->Agent::hatchnewcell(1, stem);
						//this->die();
						return;
					}
					}
			}

	/* -------------------------------------------------------------------------- */
	/*                              Differentiation                               */
	/* -------------------------------------------------------------------------- */
	#ifdef MODEL_SCAFFOLD
		in = this->index[read_t];
		if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
			if (fmod((float)hours, 48 == 0)) { // differentiation attempted every 48 hours
				Stem::differentiateStem(1, progen);
			}
		}
	#endif

	/* -------------------------------------------------------------------------- */
	/*                                  MOVEMENT                                  */
	/* -------------------------------------------------------------------------- */
   
	// Activated chondrocytes only move along their preferred gradient if there is damage.	
	//	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

	if (totaldamage != 0) this->cellSniff();

	/* -------------------------------------------------------------------------- */
	/*                      ECM PROTEIN & CHEMICAL SYNTHESIS                      */
	/* -------------------------------------------------------------------------- */
	// Calculates chemical gradients and patch chemical concentrations:
	//	int in = this->index[read_t];
	float meanTNF = this->meanNeighborChem(pTNF);
	float meanTGF = this->meanNeighborChem(pTGF);
	float meanIL1 = this->meanNeighborChem(pIL1beta);
	int countnHA = 0;
	int countfHA = 0;
	float patchTNF = this->agentWorldPtr->WHWorldChem.pTNF[in];
	float patchTGF = this->agentWorldPtr->WHWorldChem.pTGF[in];
	float patchIL1beta = (this->agentWorldPtr->WHWorldChem.pIL1beta[in]);

	// Makes collagen and aggrecan every 12 hours.
	#ifndef CALIBRATION
		if (fmod(((Agent::agentWorldPtr)->reportHour()), 1) == 0) {
	#else 
		if (fmod(((Agent::agentWorldPtr)->reportHour()), 1) == 0) { //12
	#endif

		#ifdef MODEL_SCAFFOLD
			// Active cell adhered to Ca-Alg synthesis ECM according the substrate mechanical properties:
			int in = this->index[read_t];
			if (Agent::agentPatchPtr[in].type[read_t] == CaAlg){
				for (int i = 0; i < Stem::collagenSynthRate; i++){
					this->makeOCollagen(meanTGF, meanIL1); 
				}

				for (int i = 0; i < Stem::aggrecanSynthRate; i++){
					this->makeOAggrecan(meanTNF, meanTGF, meanIL1); 
				}

			} else {
				this->makeOCollagen(meanTGF, meanIL1); 
				this->makeOAggrecan(meanTNF, meanTGF, meanIL1); 
			} 

		#else
				this->makeOCollagen(meanTGF, meanIL1); 
				this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
		#endif
		}

	// Change in chemicals due to cells:
	#ifndef CALIBRATION
		(this->agentWorldPtr->WHWorldChem.dTGF[in]) += Stem::cytokineSynthesis[0] + Cell::cytokineSynthesis[1]*(patchTGF) + Cell::cytokineSynthesis[2]*(patchIL1beta) + Cell::cytokineSynthesis[3]*(patchTNF);			//(this->agentWorldPtr->WHWorldChem.dTGF[in]) +=  Chondrocyte::cytokineSynthesis[0] + Chondrocyte::cytokineSynthesis[1]*(1 + Chondrocyte::cytokineSynthesis[2]*patchTNF);
		(this->agentWorldPtr->WHWorldChem.dTNF[in]) += Stem::cytokineSynthesis[1] + (Cell::cytokineSynthesis[5]*((patchIL1beta)/(1 + Cell::cytokineSynthesis[6]*patchTGF));	//(this->agentWorldPtr->WHWorldChem.dTNF[in]) += Chondrocyte::cytokineSynthesis[3] + Chondrocyte::cytokineSynthesis[4]/(1 + patchTGF*Chondrocyte::cytokineSynthesis[5]);
		(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += Stem::cytokineSynthesis[2] + (Cell::cytokineSynthesis[8]* ((patchTNF)/(1 + Cell::cytokineSynthesis[9]*patchTGF)); //(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += Chondrocyte::cytokineSynthesis[6] + (Chondrocyte::cytokineSynthesis[7]*patchTNF)/(Chondrocyte::cytokineSynthesis[8] + Chondrocyte::cytokineSynthesis[9]*patchTGF);
	#else
		(this->agentWorldPtr->WHWorldChem.dTGF[in]) += 50 + (2.25*patchTGF + 1.3*patchIL1beta + 5.11*patchTNF); 				//9.98 + 2.58*patchTGF + 5.11*patchTNF;				//2.11 + 3.7*patchTGF;
		(this->agentWorldPtr->WHWorldChem.dTNF[in]) += 0 + (2.42*patchIL1beta)/(1 + 4.22*patchTGF);				//5.16 + (2.42*patchIL1beta)/(1 + 4.22*patchTGF);	//2.4*patchIL1beta + 4.8/(1 + 1.27*patchTGF);		
		(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += -0.807 + (5.43*patchTNF)/(1 + 3.26*patchTGF);		//2.11 + (5.43*patchTNF)/(1 + 3.26*patchTGF);		//4;
	#endif

	/* -------------------------------------------------------------------------- */
	/*                                DEACTIVATION                                */
	/* -------------------------------------------------------------------------- */
  	// Activated chondrocytes might be deactivated once the damage is cleared:
	/*
	totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];
	
	#ifndef CALIBRATION
		if (totaldamage == 0 && rollDice(Chondrocyte::activation[4])) this->chondDeactivation();
	#else  
		if (totaldamage == 0 && rollDice(2.5)) this->chondDeactivation();
	#endif  
	*/
	/* -------------------------------------------------------------------------- */
	/*                                    DEATH                                   */
	/* -------------------------------------------------------------------------- */
	// Stem cells in Ca-Alg hydrogel have a low apoptosis rate
	#ifdef CALIBRATION
		if (rollDice(Stem::apoptosisChance)) {
			this->die();
			return;
		}
	#else
		if (rollDice(0.05)) { // from Netlogo model
			this->die();
			return;
		}
	#endif
    	// Activated chondrocytes can die naturally:
		this->life[write_t] = this->life[read_t] - 1;
		if (this->life[read_t] <= 0) {
			this->die();
		}
}

void Progen::progen_cellFunction() {
	int in = this->index[read_t];
	double hours = Agent::agentWorldPtr->reportHour();
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

	/* -------------------------------------------------------------------------- */
	/*                                PROLIFERATION                               */
	/* -------------------------------------------------------------------------- */

#ifdef MODEL_SCAFFOLD
	in = this->index[read_t];
	if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
#ifdef CALIBRATION
		if (fmod((float)hours, Progen::proliferation[0]) == 0) {
#else 
		if (fmod(hours, 24) == 0) {
#endif 
			float meanTNF = this->meanNeighborChem(TNF);
			float meanTGF = this->meanNeighborChem(TGF);
			float meanIL1 = this->meanNeighborChem(IL1beta);
			//int countfHA = this->countNeighborECM(fha);

#ifndef CALIBRATION
			float progenProlif = log10(1 + meanTNF - meanIL1 + meanTGF);
			if (rollDice(progenProlif) {
#else  
			float progenProlif = log10(1 + meanTNF - meanIL1 + meanTGF);
				if (rollDice(progenProlif)) {
#endif  
					this->Agent::hatchnewcell(1, progen);
					//this->die();
						return;
				}
			}
			} else {
#endif

#ifndef CALIBRATION
		if (fmod((float)hours, Progen::proliferation[0])) == 0 {
#else  
		if (fmod(hours, 24) == 0) {
#endif  
			float meanTNF = this->meanNeighborChem(TNF);
			float meanTGF = this->meanNeighborChem(TGF);
			float meanIL1 = this->meanNeighborChem(IL1beta);
			//int countfHA = this->countNeighborECM(fha);

#ifndef CALIBRATION
			float progen = log10(1 + meanTNF - meanIL1 + meanTGF);
			if (rollDice(progenProlif) {
#else
			float progenProlif = log10(1 + meanTNF - meanIL1 + meanTGF);
				if (rollDice(progenProlif)) {
#endif
					this->Agent::hatchnewcell(1, progen);
					//this->die();
						return;
				}
			}
			}
	/* -------------------------------------------------------------------------- */
	/*                              Differentiation                               */
	/* -------------------------------------------------------------------------- */
#ifdef MODEL_SCAFFOLD
	in = this->index[read_t];
	if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
		if (fmod((float)hours, 48 == 0)) { // differentiation attempted every 48 hours
			Progen::differentiateProgen(1, np);
		}
	}

	/* -------------------------------------------------------------------------- */
	/*                                  MOVEMENT                                  */
	/* -------------------------------------------------------------------------- */

	// Activated chondrocytes only move along their preferred gradient if there is damage.	
	//	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

	if (totaldamage != 0) this->cellSniff();

	/* -------------------------------------------------------------------------- */
	/*                      ECM PROTEIN & CHEMICAL SYNTHESIS                      */
	/* -------------------------------------------------------------------------- */
	// Calculates chemical gradients and patch chemical concentrations:
	//	int in = this->index[read_t];
	float meanTNF = this->meanNeighborChem(pTNF);
	float meanTGF = this->meanNeighborChem(pTGF);
	float meanIL1 = this->meanNeighborChem(pIL1beta);
	int countnHA = 0;
	int countfHA = 0;
	float patchTNF = this->agentWorldPtr->WHWorldChem.pTNF[in];
	float patchTGF = this->agentWorldPtr->WHWorldChem.pTGF[in];
	float patchIL1beta = (this->agentWorldPtr->WHWorldChem.pIL1beta[in]);

	// Makes  aggrecan every 12 hours.
#ifndef CALIBRATION
	if (fmod(((Agent::agentWorldPtr)->reportHour()), 1) == 0) {
#else 
	if (fmod(((Agent::agentWorldPtr)->reportHour()), 1) == 0) { //12
#endif

#ifdef MODEL_SCAFFOLD
		// Active cell adhered to Ca-Alg synthesis ECM according the substrate mechanical properties:
		int in = this->index[read_t];
		if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
			for (int i = 0; i < Progen::aggrecanSynthRate; i++) {
				this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
			}

		}
		else {
			this->makeOCollagen(meanTGF, meanIL1);
			this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
		}

#else
		this->makeOCollagen(meanTGF, meanIL1);
		this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
#endif
	}

	// Change in chemicals due to cells:
#ifndef CALIBRATION
	(this->agentWorldPtr->WHWorldChem.dTGF[in]) += Progen::cytokineSynthesis[0] + Cell::cytokineSynthesis[1]*(patchTGF) + Cell::cytokineSynthesis[2]*(patchIL1beta) + Cell::cytokineSynthesis[3]*(patchTNF);			//(this->agentWorldPtr->WHWorldChem.dTGF[in]) +=  Chondrocyte::cytokineSynthesis[0] + Chondrocyte::cytokineSynthesis[1]*(1 + Chondrocyte::cytokineSynthesis[2]*patchTNF);
	(this->agentWorldPtr->WHWorldChem.dTNF[in]) += Progen::cytokineSynthesis[1] + (Cell::cytokineSynthesis[5] * ((patchIL1beta) / (1 + Cell::cytokineSynthesis[6] * patchTGF));	//(this->agentWorldPtr->WHWorldChem.dTNF[in]) += Chondrocyte::cytokineSynthesis[3] + Chondrocyte::cytokineSynthesis[4]/(1 + patchTGF*Chondrocyte::cytokineSynthesis[5]);
	(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += Progen::cytokineSynthesis[2] + (Cell::cytokineSynthesis[8] * ((patchTNF) / (1 + Cell::cytokineSynthesis[9] * patchTGF)); //(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += Chondrocyte::cytokineSynthesis[6] + (Chondrocyte::cytokineSynthesis[7]*patchTNF)/(Chondrocyte::cytokineSynthesis[8] + Chondrocyte::cytokineSynthesis[9]*patchTGF);
#else
	(this->agentWorldPtr->WHWorldChem.dTGF[in]) += 1 + (2.25 * patchTGF + 1.3 * patchIL1beta + 5.11 * patchTNF); 				//9.98 + 2.58*patchTGF + 5.11*patchTNF;				//2.11 + 3.7*patchTGF;
	(this->agentWorldPtr->WHWorldChem.dTNF[in]) += 2.58 + (2.42 * patchIL1beta) / (1 + 4.22 * patchTGF);				//5.16 + (2.42*patchIL1beta)/(1 + 4.22*patchTGF);	//2.4*patchIL1beta + 4.8/(1 + 1.27*patchTGF);		
	(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += 0 + (5.43 * patchTNF) / (1 + 3.26 * patchTGF);		//2.11 + (5.43*patchTNF)/(1 + 3.26*patchTGF);		//4;
#endif

	/* -------------------------------------------------------------------------- */
	/*                                DEACTIVATION                                */
	/* -------------------------------------------------------------------------- */
	// Activated chondrocytes might be deactivated once the damage is cleared:
	/*
	totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

	#ifndef CALIBRATION
		if (totaldamage == 0 && rollDice(Chondrocyte::activation[4])) this->chondDeactivation();
	#else
		if (totaldamage == 0 && rollDice(2.5)) this->chondDeactivation();
	#endif
	*/
	/* -------------------------------------------------------------------------- */
	/*                                    DEATH                                   */
	/* -------------------------------------------------------------------------- */
	// Progenitor cells in Ca-Alg hydrogel have a low apoptosis rate
#ifdef CALIBRATION
	if (rollDice(Progen::apoptosisChance)) {
		this->die();
		return;
	}
#else
	if (rollDice(0.1)) { // from Netlogo model
		this->die();
		return;
	}
#endif
	// can die naturally:
	this->life[write_t] = this->life[read_t] - 1;
	if (this->life[read_t] <= 0) {
		this->die();
	}
#endif
}

void Cell::makeOCollagen(float meanTGF, float meanIL1) {
	int read_index;

	// Check if the location has been modified in this tick
	if (isModified(this->index)) {
		read_index = write_t;  // If it has, work off of the intermediate value
	}
	else {
		read_index = read_t;  // If it has NOT, work off of the original value
	}

	int dx, dy, dz;

	// Location of cell in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];

	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;
	int randInt, target, in;
	vector <int> neighbors;
	//vector <int> damagedneighbors;

	// Make a list of neighboring patches
	#ifndef MODEL_3D
	for (int i = 9; i < 18; i++) {
	#else
	for (int i = 0; i < 27; i++) {
	#endif
		dx = Agent::dX[i];
		dy = Agent::dY[i];
		dz = Agent::dZ[i];
		in = (x + dx) + (y + dy) * nx + (z + dz) * nx * ny;

		// Try a new neighboring patch if this one is outside the world dimensions:
		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;

		// Add the valid neighboring patch to the list:
		neighbors.push_back(i);
	}

	// Target a random damaged neighboring patch, if there are any.
	if (neighbors.size() > 0) {
		int tid = 0;
	#ifdef _OMP
		tid = omp_get_thread_num();     // Get thread id in order to access the seed that belongs to this thread
	#endif

		randInt = rand_r(&(agentWorldPtr->seeds[tid])) % neighbors.size();
		target = neighbors[randInt];
		dx = Agent::dX[target];
		dy = Agent::dY[target];
		dz = Agent::dZ[target];

		// Move to new patch and sprout ocollagen
		in = (x + dx) + (y + dy) * nx + (z + dz) * nx * ny;
		this->move(dx, dy, dz, read_index);

		Agent::agentECMPtr[in].ocollagen[write_t] = Agent::agentECMPtr[in].ocollagen[read_t] + 1 + rand() % 2;
	#ifdef OPT_ECM
		Agent::agentECMPtr[in].set_dirty();
	#endif
	}

	//// Make a list of damaged neighboring patches
	//#ifndef MODEL_3D
	//	for (int i = 9; i < 18; i++) {
	//#else
	//	for (int i = 0; i < 27; i++) {
	//#endif
	//		dx = Agent::dX[i];
	//		dy = Agent::dY[i];
	//		dz = Agent::dZ[i];
	//		in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
 //   
	//		// Try a new neighboring patch if this one is outside the world dimensions:
	//		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
	//		
	//		// Add the valid damaged neighboring patch to the list:
	//		if (Agent::agentPatchPtr[in].damage[read_t] != 0) damagedneighbors.push_back(i);
	//	}

	//// Target a random damaged neighboring patch, if there are any.
	//if (damagedneighbors.size() > 0) {
	//	int tid = 0;
	//	#ifdef _OMP
	//		tid = omp_get_thread_num();     // Get thread id in order to access the seed that belongs to this thread
	//	#endif

	//	randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();
	//	target = damagedneighbors[randInt];
	//	dx = Agent::dX[target];
	//	dy = Agent::dY[target];
	//	dz = Agent::dZ[target];

	//	// Based on chance, chemical concentrations, move to new patch and sprout ocollagen
	//	#ifndef CALIBRATION
	//		int stimulation = (Cell::ECMsynthesis[0]*(log10(meanTGF + Cell::ECMsynthesis[1])/(Cell::ECMsynthesis[2] + meanIL1*Cell::ECMsynthesis[2]));
	//		if ((rollDice(Cell::ECMsynthesis[3] + stimulation)) || (rollDice(Cell::ECMsynthesis[4] + stimulation/Cell::ECMsynthesis[5])) || (rollDice(Cell::ECMsynthesis[6] + stimulation/Cell::ECMsynthesis[7]))) {
	//	#else 
	//		int stimulation = ((log10(1 + (meanTGF)))/(1+ meanIL1)); 
	//		if ((rollDice(50 + stimulation)) || (rollDice(25 + stimulation/2)) || (rollDice(10+stimulation/5))) {
	//	#endif 
	//			in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
	//			this->move(dx, dy, dz, read_index);
	//			
	//			Agent::agentECMPtr[in].ocollagen[write_t] = Agent::agentECMPtr[in].ocollagen[read_t] + 1 + rand()%2;
	//			#ifdef OPT_ECM
	//				Agent::agentECMPtr[in].set_dirty(); 
	//			#endif
	//		}
	//	}
}

void Cell::makeOAggrecan(float meanTNF, float meanTGF, float meanIL1) {
	int read_index;
	// Check if the location has been modified in this tick:
	if (isModified(this->index)) {
		read_index = write_t;		// If it has, work off of the intermediate value
	} else {
		read_index = read_t;		// If it has NOT, work off of the original value
	}

	int dx, dy, dz;

  	// Location of cell in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];

  	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;
	int randInt, target, in;
	vector <int> neighbors;

	// Make a list of neighboring patches
	#ifndef MODEL_3D
	for (int i = 9; i < 18; i++) {
	#else
	for (int i = 0; i < 27; i++) {
	#endif
		dx = Agent::dX[i];
		dy = Agent::dY[i];
		dz = Agent::dZ[i];
		in = (x + dx) + (y + dy) * nx + (z + dz) * nx * ny;

		// Try a new neighboring patch if this one is outside the world dimensions:
		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;

		// Add the valid neighboring patch to the list:
		neighbors.push_back(i);
	}
	
	// Target a random neighboring patch, if there are any.
	if (neighbors.size() > 0) {
		int tid = 0;
#ifdef _OMP
		tid = omp_get_thread_num();     // Get thread id in order to access the seed that belongs to this thread
#endif

		randInt = rand_r(&(agentWorldPtr->seeds[tid])) % neighbors.size();
		target = neighbors[randInt];
		dx = Agent::dX[target];
		dy = Agent::dY[target];
		dz = Agent::dZ[target];

		// Move to new patch and sprout ocollagen
		in = (x + dx) + (y + dy) * nx + (z + dz) * nx * ny;
		this->move(dx, dy, dz, read_index);

		Agent::agentECMPtr[in].ocollagen[write_t] = Agent::agentECMPtr[in].ocollagen[read_t] + 1 + rand() % 2;
#ifdef OPT_ECM
		Agent::agentECMPtr[in].set_dirty();
#endif
	}
}

//void Chondrocyte::makeHyaluronan(float meanTNF, float meanTGF, float meanIL1) {
//	int read_index;
//	// Check if the location has been modified in this tick
//	if (isModified(this->index)) {
//		read_index = write_t;		// If it has, work off of the intermediate value
//	} else {
//		read_index = read_t;		// If it has NOT, work off of the original value
//	}
//	
//  	// Location of chondrocyte in x,y,z dimensions of world.
//	int x = this->ix[read_index];
//	int y = this->iy[read_index];
//	int z = this->iz[read_index];
//  
//  	// Number of patches in x,y,z dimensions of world
//	int nx = Agent::nx;
//	int ny = Agent::ny;
//	int nz = Agent::nz;
//	vector <int> damagedneighbors;
//
//	#ifndef CALIBRATION
//		int stimulation = Chondrocyte::ECMsynthesis[12]*(log10(Chondrocyte::ECMsynthesis[12] + meanTGF + meanTNF + meanIL1)) + Chondrocyte::ECMsynthesis[13];
//	#else  
//		int stimulation = 1*log10(1 + meanTGF + meanTNF + meanIL1) + 0;
//	#endif  
//
//	int dx, dy, dz, in, randInt, target;
//
//	// Make list of damaged neighboring patches
//	#ifndef MODEL_3D
//		for (int i = 9; i < 18; i++) {
//	#else
//		for (int i = 0; i < 27; i++) {
//	#endif
//			dx = Agent::dX[i];
//			dy = Agent::dY[i];
//			dz = Agent::dZ[i];
//			in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
//			
//			// Try a new neighboring patch if this one is outside the world dimensions.
//			if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
//		
//			// Add the valid damaged neighboring patch to the list.
//			if (Agent::agentPatchPtr[in].damage[read_t] != 0) damagedneighbors.push_back(i);
//	}
//
//	// Target a random damaged neighboring patch, if there are any.
//	if (damagedneighbors.size() > 0) {
//		int tid = 0;
//		#ifdef _OMP
//			tid = omp_get_thread_num();    // Get thread id in order to access the seed that belongs to this thread
//		#endif
//
//		randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();
//		target = damagedneighbors[randInt];
//		dx = Agent::dX[target];
//		dy = Agent::dY[target];
//		dz = Agent::dZ[target];
//		in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
//
//		// Based on mean TGF, TNF, IL1, chance, move to new patch and sprout HA
//		#ifndef CALIBRATION
//			if (rollDice(Chondrocyte::ECMsynthesis[14] + stimulation)) { 			
//		#else  
//			if (rollDice(50 + stimulation)) { 
//		#endif  
//				this->move(dx, dy, dz, read_index);
//
//			}
//		#ifndef CALIBRATION
//			} else if (rollDice(Chondrocyte::ECMsynthesis[15] + stimulation/Chondrocyte::ECMsynthesis[16])) {
//		#else  
//			} else if (rollDice(5 + stimulation/10)) {		
//		#endif 
//				in = this->index[read_t];
//				#ifdef OPT_ECM
//					Agent::agentECMPtr[in].set_dirty(); 
//				#endif
//			}
//}

void Stem::differentiateStem(int number = 1, int agentType = progen) {
	// stem cells differentiate to the next stage, progenitor

	float meanTNF = this->meanNeighborChem(TNF);
	float meanTGF = this->meanNeighborChem(TGF);
	float meanIL1 = this->meanNeighborChem(IL1beta);

	int in = this->index[read_t];
	float stemDiff = 0.5 + (Stem::differentiation[3] * meanTGF);
	if (rollDice(stemDiff)) {
		if (rollDice(0.7)) { // asymmetric differentiation
			Agent::hatchnewcell(number, agentType); // only create new progenitor cell nearby
		}
		else { // symmetric differentiation
			Cell* temp = this;
			Progen* dp2 = dynamic_cast<Progen*>(temp); // dynamic cast to next derived class (cell type)
			if (dp2 == nullptr) {
				cout << "Casting Failed" << endl;
			}
			Agent::agentPatchPtr[in].occupiedby[write_t] = agentType; // switch current patch to be marked as occupied by new cell type
			Agent::hatchnewcell(number, agentType); // create new progenitor cell nearby
		}
	}
}

void Progen::differentiateProgen(int number = 1, int agentType = np) {
	// np progenitor cells differentiate to the next stage, np cells

	int in = this->index[read_t];
	if (rollDice(0.7)) { // asymmetric differentiation
		Agent::hatchnewcell(number, agentType); // only create new NP cell nearby
	}
	else { // symmetric differentiation
		Cell* temp = this;
		NP* dp2 = dynamic_cast<NP*>(temp); // dynamic cast to next derived class (cell type)
		if (dp2 == nullptr) {
			cout << "Casting Failed" << endl;
		}
		else {
			cout << "Casting Successful" << endl;
		}
		Agent::agentPatchPtr[in].occupiedby[write_t] = agentType; // switch current patch to be marked as occupied by new cell type
		Agent::hatchnewcell(number, agentType); // create new NP cell nearby
	}
}
