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
float Cell::cytokineSynthesis[10] = {10, 0.05, 10, 5, 2.4, 4, 2, 5, 1, 3.2}; //calibration variables
float Cell::activation[5] = {500, 50, 0, 25, 2.5};
float Cell::ECMsynthesis[12] = {1, 1, 1, 50, 25, 2, 10, 5, 1, 1, 25, 2};
float Cell::proliferation[6] = {24, 10, 1, 0, 25, 2};

int Stem::numOfStem = 0;
float Stem::migrationSpeed = 1; // patch/tick

float Stem::CaAlgMigration[2] = { 0.11, 0.35 };
float Stem::cytokineSynthesis[3] = { 50, 0, -0.807 };
float Stem::CollagenSynth[3] = {  };
float Stem::AggrecanSynth[3] = {  };
float Stem::ECMsynthesis[4] = {};
float Stem::proliferation[4] = {};
float Stem::differentiation[5] = {};

int Progen::numOfProgen = 0; 
float Progen::migrationSpeed = 1;    // patch/tick

float Progen::CaAlgMigration[2] = { 0.11, 0.83 };
float Progen::cytokineSynthesis[3] = {  };
float Progen::ECMsynthesis[4] = {};
float Progen::proliferation[4] = {};
float Progen::differentiation[5] = {};

int NP::numOfNP = 0;
float NP::migrationSpeed = 1;    // patch/tick
float collagenSynthRate = 1;
float aggrecanSynthRate = 1.5;

float NP::CaAlgMigration[2] = { 0.11, 1.30 };
float Stem::CollagenSynth[3] = { 10, 6.45, 3.6 };
float Stem::AggrecanSynth[3] = { 10, 38, 16.6 };



Cell::Cell() {
	cout << "default cell alloc" << endl;
}

Cell::Cell(Patch* patchPtr) {
	this->ix[write_t] = patchPtr->indice[0];
	this->iy[write_t] = patchPtr->indice[1];
	this->iz[write_t] = patchPtr->indice[2];
	this->index[write_t] = patchPtr->index;
	this->alive[write_t] = true;

	#ifndef MODEL_SCAFFOLD
		// Unactivated chondrocytes live for 5 to 11 days. 0 corresponds to hours:
		if (Agent::agentWorldPtr->clock == 0) this->life[write_t] = WHWorld::reportTick(0, rand()%12);
		else this->life[write_t] = WHWorld::reportTick(0, 5 + rand()%7);
	#else
		this->life[write_t] = WHWorld::reportTick(0, 5 + rand()%7);
	#endif

	this->activate[write_t] = false;
	this->color[write_t] = ccell;
	this->size[write_t] = 2;
	this->type[write_t] = cell;
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
	this->color[read_t] = ccell;
	this->size[read_t] = 2;
	this->type[read_t] = cell;
	/* Added by MM to check types of cell stages and add to respective counters: */
	if (typeid(*this) == typeid(Stem)) {
		Stem::numOfStem++;
	} else if (typeid(*this) == typeid(Progen)) {
		Progen::numOfProgen++:
	} else if (typeid(*this) == typeid(NP)) {
		NP::numOfNP++;
	}
	Cell::numOfCells++;
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
		Progen::numOfProgen++:
	}
	else if (typeid(*this) == typeid(NP)) {
		NP::numOfNP++;
	}
	Cell::numOfCells++;  
}

Cell::~Cell() {}

void Cell::cellFunction() {
	if (this->alive[read_t] == false) return;
	//if (this->activate[read_t] == false) this->chond_cellFunction();
	//else this->achond_cellFunction();
}

void Cell::chond_cellFunction() {  //NOTE MM: this function code has to be copied into stem_cellFunction, progen_cellFunction, and NP_cellFunction. 
								   //chond_cellFunction will no longer be used
	int in = this->index[read_t];
	double hours = Agent::agentWorldPtr->reportHour();
	int totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];

	/* Unactivated chondrocytes only move along their preferred gradient and perform biological functions if there is damage. */
	if (totaldamage == 0) {

		#ifdef MODEL_SCAFFOLD
			// If cell is moving on Ca-Alg substrate chondrocyte move up to "migrationSpeed" patches per tick determined by substrate composition
			if (Agent::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
				for (int dx = 0; dx < Agent::migrationSpeed; dx++) this->wiggle(); 
			
			} else if(rollDice(0.25)){ // If cell is not actively migrating, consider chance of moving to next patch            
				this->wiggle(); 
			}
		#else
				this->wiggle();
		#endif

	} else {

		/* -------------------------------------------------------------------------- */
		/*                                PROLIFERATION                               */
		/* -------------------------------------------------------------------------- */
		#ifdef MODEL_SCAFFOLD
			// Cells in CaAlg hydrogel proliferate at a proliferation rate given time and hydrogel composition (% Alg) 
			in = this->index[read_t];
			if (Agent::agentPatchPtr[in].type[read_t] == CaAlg){

			#ifdef CALIBRATION
				if ((fmod((float)Agent::agentWorldPtr->clock, 2) == 0) && rollDice (Agent::proliferationRate)){   //check every hour 
			#else 
				if ((fmod((float)Agent::agentWorldPtr->clock, 2) == 0) && rollDice (Agent::proliferationRate)){ //check every hour
			#endif

				float meanTNF = this->meanNeighborChem(TNF);
				float meanTGF = this->meanNeighborChem(TGF);
				float meanIL1 = this->meanNeighborChem(IL1beta);
				int countfHA = 0;
				int TGFrelated = 0;

			#ifndef CALIBRATION
				if (meanTGF <= 10) { 
			#else 
				if (meanTGF <= 10) {
			#endif
					TGFrelated = 1;  // Low TGF (0.1-1 ng) stimulate chond proliferation and attraction. 
				} else {
					TGFrelated = -1; // High TGF (1-10 ng) inhibits proliferation.
				}

			#ifndef CALIBRATION
				float chondProlif = Chondrocyte::proliferation[2]*(log10(1 + meanTNF + meanIL1 + TGFrelated*meanTGF)) + Chondrocyte::proliferation[3];
				if (rollDice(Chondrocyte::proliferation[4] + chondProlif/Chondrocyte::proliferation[5])) {
			#else
				float chondProlif = 1*log10(1 + meanTNF + meanIL1 + TGFrelated*meanTGF) + 0;
				if (rollDice(25 + chondProlif/2)) { 
			#endif	
					this->hatchnewchondrocyte(2);
					this->die();
					return;
				}
				}				

				// Unactivated chondrocytes proliferate every 24 hours.
				} else{
		#endif

	#ifndef CALIBRATION
		if (fmod((float) hours, Chondrocyte::proliferation[0]) == 0) {
	#else  
		if (fmod(hours, 24) == 0) {
	#endif  
			float meanTNF = this->meanNeighborChem(TNF);
			float meanTGF = this->meanNeighborChem(TGF);
			float meanIL1 = this->meanNeighborChem(IL1beta);
			int countfHA =  this->countNeighborECM(fha);
			int TGFrelated = 0;

		#ifndef CALIBRATION
			if (meanTGF <= 10) {
		#else
			if (meanTGF <= 10) {
		#endif
				TGFrelated = 1;  // Low TGF (0.1-1 ng) stimulate chond proliferation and attraction.
			} else {
				TGFrelated = -1;  // High TGF (1-10 ng) inhibits proliferation.
			}

		#ifndef CALIBRATION
			float chondProlif = Chondrocyte::proliferation[2]*(log10(1 + meanTNF + meanIL1 + TGFrelated*meanTGF)) + Chondrocyte::proliferation[3];  
			if (rollDice(Chondrocyte::proliferation[4] + chondProlif/Chondrocyte::proliferation[5])) {  
		#else  
			float chondProlif = log10(1 + meanTNF + meanIL1 + TGFrelated*meanTGF);
			if (rollDice(25 + chondProlif/2)) {		
		#endif 
				this->hatchnewchondrocyte(2);
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
	// An unactivated chondrocyte can be activated if it is in the damage zone:
	if (Agent::agentPatchPtr[in].inDamzone == true) {
		// Low TGF promote and high TGF inhibit chances of chondrocyte activation:
		int patchTGF = agentWorldPtr->WHWorldChem.pTGF[in];
		
	#ifndef CALIBRATION
		if ((patchTGF > 10 && rollDice(Chondrocyte::activation[1])) || (patchTGF > Chondrocyte::activation[2]) || (rollDice(Chondrocyte::activation[3]))) { 
	#else  
		if ((patchTGF > 10 && rollDice(50.0)) || (patchTGF > 0) || (rollDice(25))) { 
	#endif
			this->chondActivation();
		}
		}
	}

	/* -------------------------------------------------------------------------- */
	/*                                    DEATH                                   */
	/* -------------------------------------------------------------------------- */
	
	// Cells in Ca-Alg hydrogel have at a viability/death rate determined by time:
	#ifdef MODEL_SCAFFOLD
		#ifdef CALIBRATION
			if (fmod((float)Agent::agentWorldPtr->clock, Agent::CaAlgViability[2]) == 0 && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
		#else
			if (fmod((float)Agent::agentWorldPtr->clock, 6) && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
		#endif
				if (rollDice(100 - Agent::viabilityRate)){
					this->die();
					return; 
				}
			}
	#endif

    // Unactivated chondrocytes can die naturally:
	this->life[write_t] = this->life[read_t] - 1;
	if (this->life[read_t] <= 0) this->die();	
}

void Cell::hatchnewcell(int number, int agentType) {
	int newcells = 0;
  
  	// Location of cell in x,y,z dimensions of the world
	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];
  
  	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	// Shuffle neighboring patches and go through them in a random order:
	#ifdef MODEL_3D
		random_shuffle(&Agent::neighbor[0], &Agent::neighbor[26]);
		for (int i = 0; i < 27 && newcells < number; i++) {
	#else
		random_shuffle(&Agent::neighbor[0], &Agent::neighbor[7]);
		for (int i = 0; i < 8 && newcells < number; i++) {
	#endif

		// Distance away from target neighboring patch in x,y,z dimensions
		int dx = Agent::dX[neighbor[i]];
		int dy = Agent::dY[neighbor[i]];
		int dz = Agent::dZ[neighbor[i]];
			
		// Patch row major index of target neighboring patch:
		int in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

		// Try a new target neighboring patch if this one is not inside the world dimensions, or is occupied.
		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
		int targetType = agentPatchPtr[in].type[read_t];

		// Create a new cell of agentType at the valid target neighboring patch:
		Cell* newcell = NULL;
		switch (agentType) {
			case stem:
				Stem* newcell = new Stem(x + dx, y + dy, z + dz);
				break;
			case progen:
				Progen* newcell = new Progen(x + dx, y + dy, z + dz);
				break;
			case np:
				NP* newcell = new NP(x + dx, y + dy, z + dz);
				break;
		}
		newcells++;

		// Update target neighboring patch as occupied:
		Agent::agentPatchPtr[in].setOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = agentType;
			
		/* If executing OMP version, add the pointer to this new cell to the thread-local list first. 
			* WHWorld::UpdateCells() will take care of putting it in the global list at the end 
			*/
		#ifdef _OMP
			int tid = omp_get_thread_num();
			Agent::agentWorldPtr->localNewCells[tid]->push_back(newcell);
		#else
			// If executing serial version, add the pointer to this new cell to the global list right away
			Agent::agentWorldPtr->cells.addData(newcell, DEFAULT_TID);
		#endif

		newcell->wiggle();
	}
}

void Cell::cellSniff() {
	int in = this->index[read_t];
	if ((Agent::agentPatchPtr[in]).inDamzone == true) {
		if (rollDice(80) && this->moveToHighestChem(pcellgrad) == true){

			#ifdef MODEL_SCAFFOLD
				// If cell is moving on Ca-Alg substrate:
				if (Agent::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
					
					// Move up to "migrationSpeed" patches per tick:
					for (int dx = 0; dx < Agent::migrationSpeed; dx++) this->wiggle(); 

				// If cell is not actively migrating, consider chance of moving to next patch:
				} else if(rollDice(0.25)){		
					this->wiggle(); 
				}
			#else
					this->wiggle();
			#endif
        }

	} else {
		#ifdef MODEL_SCAFFOLD
			// If cell is moving on Ca-Alg substrate:
			if (Agent::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
				for (int dx = 0; dx < Agent::migrationSpeed; dx++) this->wiggle(); // Move up to "migrationSpeed" patches per tick:
			} else if(rollDice(0.25)) this->wiggle(); 							   // If cell is not actively migrating, consider chance of moving to next patch
		#else
			this->wiggle();
		#endif
	}

	// TGF can excite stem cell and overcome gradient: //NOTE MM: double check this for MSCs
	if ((this->meanNeighborChem(pTGF)) > 0) {
		#ifdef MODEL_SCAFFOLD
			// If cell is moving on Ca-Alg substrate:
			if (Agent::migrationSpeed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
				for (int dx = 0; dx < Agent::migrationSpeed; dx++) this->wiggle(); // Move up to "migrationSpeed" patches per tick:
			} else if(rollDice(0.25)) this->wiggle(); 							   // If cell is not actively migrating, consider chance of moving to next patch: 		
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

//void Chondrocyte::chondDeactivation() {
//	int in = this->index[read_t];
//	int target = this->index[write_t]; // This assumes that after this function is called, no more move() would be called in the same tick 
//	if (this->activate[read_t] == true&& this->life[read_t] > 1) {
//		this->type[write_t] = chondrocyte;
//		this->activate[write_t] = false;
//		this->color[write_t] = cchondrocyte;
//		Agent::agentPatchPtr[in].setOccupied();
//		Agent::agentPatchPtr[in].occupiedby[write_t] = chondrocyte;
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

						TGFrelated = 1;  // Low TGF (0.1-1nm) stimulate proliferation and attraction. 
					} else {
						TGFrelated = -1; // High TGF (1-10nm) inhibits proliferation. 
					}

					#ifndef CALIBRATION
						float stemProlif = log10(1 + Stem::proliferation[2]*meanTNF + Stem::proliferation[3]*meanIL1 + TGFrelated*meanTGF);
						if (rollDice(stemProlif) {  
					#else  
						float stemProlif = log10(1 + meanTNF + meanIL1 + TGFrelated*meanTGF); 
						if (rollDice(stemProlif)) { 										
					#endif  
							this->hatchnewcell(1);
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
			if (meanTGF <= meanTGF <= Stem::proliferation[0]) {
  		#else  
			if (meanTGF <= 10) {
  		#endif  
				TGFrelated = 1; // Low TGF (0.1-1nm) stimulate chond proliferation and attraction. 
			} else {
				TGFrelated = -1;  // High TGF (1-10nm) inhibits proliferation. 
			}

		#ifndef CALIBRATION
			float stemProlif = log10(1 + Stem::proliferation[2] * meanTNF + Stem::proliferation[3] * meanIL1 + TGFrelated * meanTGF);
			if (rollDice(stemProlif) {
		#else
			float stemProlif = log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF);
				if (rollDice(stemProlif)) {
		#endif
				this->hatchnewcell(1);
				//this->die();
				return;
			}
			}

	#ifdef MODEL_SCAFFOLD
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
				for (int i = 0; i < Agent::collagenSynthRate; i++){
					this->makeOCollagen(meanTGF, meanIL1); 
				}

				for (int i = 0; i < Agent::aggrecanSynthRate; i++){
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
		(this->agentWorldPtr->WHWorldChem.dTGF[in]) += Chondrocyte::cytokineSynthesis[0] + Chondrocyte::cytokineSynthesis[1]*(patchTGF + Chondrocyte::cytokineSynthesis[2]*patchTNF);			//(this->agentWorldPtr->WHWorldChem.dTGF[in]) +=  Chondrocyte::cytokineSynthesis[0] + Chondrocyte::cytokineSynthesis[1]*(1 + Chondrocyte::cytokineSynthesis[2]*patchTNF);
		(this->agentWorldPtr->WHWorldChem.dTNF[in]) += Chondrocyte::cytokineSynthesis[3] + (Chondrocyte::cytokineSynthesis[4]*patchIL1beta)/(1 + Chondrocyte::cytokineSynthesis[5]*patchTGF);	//(this->agentWorldPtr->WHWorldChem.dTNF[in]) += Chondrocyte::cytokineSynthesis[3] + Chondrocyte::cytokineSynthesis[4]/(1 + patchTGF*Chondrocyte::cytokineSynthesis[5]);
		(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += Chondrocyte::cytokineSynthesis[6] + (Chondrocyte::cytokineSynthesis[7]*patchTNF)/(Chondrocyte::cytokineSynthesis[8] + Chondrocyte::cytokineSynthesis[9]*patchTGF); //(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += Chondrocyte::cytokineSynthesis[6] + (Chondrocyte::cytokineSynthesis[7]*patchTNF)/(Chondrocyte::cytokineSynthesis[8] + Chondrocyte::cytokineSynthesis[9]*patchTGF);
	#else
		(this->agentWorldPtr->WHWorldChem.dTGF[in]) += 10 + 0.05*(patchTGF + 10*patchTNF); 				//9.98 + 2.58*patchTGF + 5.11*patchTNF;				//2.11 + 3.7*patchTGF;
		(this->agentWorldPtr->WHWorldChem.dTNF[in]) += 5 + (2.4*patchIL1beta)/(1 + 4*patchTGF);				//5.16 + (2.42*patchIL1beta)/(1 + 4.22*patchTGF);	//2.4*patchIL1beta + 4.8/(1 + 1.27*patchTGF);		
		(this->agentWorldPtr->WHWorldChem.dIL1beta[in]) += 2 + (5*patchTNF)/(1 + 3.2*patchTGF);		//2.11 + (5.43*patchTNF)/(1 + 3.26*patchTGF);		//4;
	#endif

	/* -------------------------------------------------------------------------- */
	/*                                DEACTIVATION                                */
	/* -------------------------------------------------------------------------- */
  	// Activated chondrocytes might be deactivated once the damage is cleared:
	totaldamage = ((Agent::agentWorldPtr)->worldPatch)->numOfEachTypes[damage];
	
	#ifndef CALIBRATION
		if (totaldamage == 0 && rollDice(Chondrocyte::activation[4])) this->chondDeactivation();
	#else  
		if (totaldamage == 0 && rollDice(2.5)) this->chondDeactivation();
	#endif  

	/* -------------------------------------------------------------------------- */
	/*                                    DEATH                                   */
	/* -------------------------------------------------------------------------- */
	// Cells in Ca-Alg hydrogel have at a viability/death rate determined by time 
	#ifdef CALIBRATION
		if (fmod((float)Agent::agentWorldPtr->clock, Agent::CaAlgViability[2]) == 0 && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
	#else
		if (fmod((float)Agent::agentWorldPtr->clock, 6) && Agent::agentPatchPtr[in].type[read_t] == CaAlg){
	#endif
			if (rollDice(100 - Agent::viabilityRate)){
				this->die();
				return; 
			}
		}

    	// Activated chondrocytes can die naturally:
		this->life[write_t] = this->life[read_t] - 1;
		if (this->life[read_t] <= 0) {
			this->die();
		}
	}

void Cell::makeOCollagen(float meanTGF, float meanIL1) {
	int read_index;

	// Check if the location has been modified in this tick
	if (isModified(this->index)) {
		read_index = write_t;  // If it has, work off of the intermediate value
	} else {
		read_index = read_t;  // If it has NOT, work off of the original value
	}

	int dx, dy, dz;

  	// Location of chondrocyte in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];

  	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;
	int randInt, target, in;
	vector <int> damagedneighbors;

	// Make a list of damaged neighboring patches
	#ifndef MODEL_3D
		for (int i = 9; i < 18; i++) {
	#else
		for (int i = 0; i < 27; i++) {
	#endif
			dx = Agent::dX[i];
			dy = Agent::dY[i];
			dz = Agent::dZ[i];
			in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
    
			// Try a new neighboring patch if this one is outside the world dimensions:
			if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
			
			// Add the valid damaged neighboring patch to the list:
			if (Agent::agentPatchPtr[in].damage[read_t] != 0) damagedneighbors.push_back(i);
		}

	// Target a random damaged neighboring patch, if there are any.
	if (damagedneighbors.size() > 0) {
		int tid = 0;
		#ifdef _OMP
			tid = omp_get_thread_num();     // Get thread id in order to access the seed that belongs to this thread
		#endif

		randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();
		target = damagedneighbors[randInt];
		dx = Agent::dX[target];
		dy = Agent::dY[target];
		dz = Agent::dZ[target];

		// Based on chance, chemical concentrations, move to new patch and sprout ocollagen
		#ifndef CALIBRATION
			int stimulation = (Chondrocyte::ECMsynthesis[0]*(log10(meanTGF + Chondrocyte::ECMsynthesis[1])/(Chondrocyte::ECMsynthesis[2] + meanIL1*Chondrocyte::ECMsynthesis[2]));
			if ((rollDice(Chondrocyte::ECMsynthesis[3] + stimulation)) || (rollDice(Chondrocyte::ECMsynthesis[4] + stimulation/Chondrocyte::ECMsynthesis[5])) || (rollDice(Chondrocyte::ECMsynthesis[6] + stimulation/Chondrocyte::ECMsynthesis[7]))) {
		#else 
			int stimulation = ((log10(1 + (meanTGF)))/(1+ meanIL1)); 
			if ((rollDice(50 + stimulation)) || (rollDice(25 + stimulation/2)) || (rollDice(10+stimulation/5))) {
		#endif 
				in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
				this->move(dx, dy, dz, read_index);
				
				Agent::agentECMPtr[in].ocollagen[write_t] = Agent::agentECMPtr[in].ocollagen[read_t] + 1 + rand()%2;
				#ifdef OPT_ECM
					Agent::agentECMPtr[in].set_dirty(); 
				#endif
			}
			}
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

  	// Location of chondrocyte in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];

  	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;
	int randInt, target, in;
	vector <int> damagedneighbors;

	// Make list of damaged neighboring patches:
	#ifndef MODEL_3D
		for (int i = 9; i < 18; i++) {
	#else
		for (int i = 0; i < 27; i++) {
	#endif
			dx = Agent::dX[i];
			dy = Agent::dY[i];
			dz = Agent::dZ[i];
			in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;
			
			// Try a new neighboring patch if this one is outside the world dimensions.
			if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
			
			// Add the valid damaged neighboring patch to the list.
			if (Agent::agentPatchPtr[in].damage[read_t] != 0) damagedneighbors.push_back(i);
		}

	// Target a random damaged neighboring patch, if there are any.
	if (damagedneighbors.size() > 0) {
		int tid = 0;
		#ifdef _OMP
			tid = omp_get_thread_num(); // Get thread id in order to access the seed that belongs to this thread
		#endif

		randInt = rand_r(&(agentWorldPtr->seeds[tid])) % damagedneighbors.size();
		target = damagedneighbors[randInt];
		dx = Agent::dX[target];
		dy = Agent::dY[target];
		dz = Agent::dZ[target];
		in = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

		// Based on mean TGf, IL1, TNF, chance, move to new patch and sprout oaggrecan.
		#ifndef CALIBRATION
			if 
			int stimulation = Chondrocyte::ECMsynthesis[8]*log10((meanTGF + Chondrocyte::ECMsynthesis[9])/(Chondrocyte::ECMsynthesis[9] + meanTNF));
			if (rollDice(Chondrocyte::ECMsynthesis[10] + stimulation/Chondrocyte::ECMsynthesis[11])) {
		#else 
			int stimulation = log10((1 + meanTGF)/(1  + meanTNF));
			if (rollDice(25 + stimulation/2)){ 
		#endif
				this->move(dx, dy, dz, read_index);
				Agent::agentECMPtr[in].oaggrecan[write_t] = Agent::agentECMPtr[in].oaggrecan[read_t] + 1 + rand()%2;
				#ifdef OPT_ECM
					Agent::agentECMPtr[in].set_dirty(); 
				#endif
				// cout << " sprout new aggrecan at " << in << endl; 
			}
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

void stem_cellFunction() {



}


void progen_cellFunction() {



}

void NP_cellFunction() {



}

void Stem::differentiateStem(int number = 1, int agentType = progen) {
	// stem cells differentiate to the next stage, progenitor

	int in = this->index[read_t]
	if (rollDice(0.7)) { // asymmetric differentiation
		Cell::hatchnewcell(number, agentType); // only create new progenitor cell nearby
	} else { // symmetric differentiation
		Cell* temp = &this;
		Progen* dp2 = dynamic_cast<Progen*>(temp); // dynamic cast to next derived class (cell type)
		if (dp2 == nullptr)
			cout << "Casting Failed" << endl;
		else
			cout << "Casting Successful" << endl;
		Agent::agentPatchPtr[in].occupiedby[write_t] = agentType; // switch current patch to be marked as occupied by new cell type
		Cell::hatchnewcell(number, agentType); // create new progenitor cell nearby
}

void Progen::differentiateProgen(int number = 1, int agentType = np) {
	// np progenitor cells differentiate to the next stage, np cells

	int in = this->index[read_t]
	if (rollDice(0.7)) { // asymmetric differentiation
		Cell::hatchnewcell(number, agentType); // only create new NP cell nearby
	}
	else { // symmetric differentiation
		Cell* temp = &this;
		NP* dp2 = dynamic_cast<NP*>(temp); // dynamic cast to next derived class (cell type)
		if (dp2 == nullptr)
			cout << "Casting Failed" << endl;
		else
			cout << "Casting Successful" << endl;
		Agent::agentPatchPtr[in].occupiedby[write_t] = agentType; // switch current patch to be marked as occupied by new cell type
		Cell::hatchnewcell(number, agentType); // create new NP cell nearby
}
