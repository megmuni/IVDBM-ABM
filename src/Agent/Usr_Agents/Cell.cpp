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
float Cell::proliferation[4] = { 24, 10, 1, 0 };
float Cell::cytokineSynthesis[10] = {10, 0.05, 10, 5, 2.4, 4, 2, 5, 1, 3.2}; //calibration variables

int Stem::numOfStem = 0;
float Stem::migrationSpeed = 1; // patch/tick. Not an input parameter
float Stem::OCR = 104.65 / 2; // fmol/h/cell divided by 2 for /tick
float Stem::apoptosisChance = 0.1;
float Stem::collagenSynthRate = 1; // placeholder values which will be recalculated
float Stem::aggrecanSynthRate = 0.5; // placeholder values which will be recalculated

float Stem::CaAlgMigration[2] = { 0.11, 0.35 };
float Stem::cytokineSynthesis[3] = { 5, 0, 0 };
//float Stem::ECMsynthesis[4] = {};
float Stem::CollagenSynth[1] = { 10 };
float Stem::AggrecanSynth[1] = { 100000 };
float Stem::proliferation[4] = {10, 0.8, 0.001, 0.5};
float Stem::differentiation[5] = { 0.7, 0.5, 0.001, 48 };

int Progen::numOfProgen = 0; 
float Progen::migrationSpeed = 1;    // patch/tick. Not an input parameter
float Progen::OCR = 30.46 / 2; // fmol/h/cell divided by 2 for /tick
float Progen::apoptosisChance = 0.1;
float Progen::aggrecanSynthRate = 1;

float Progen::CaAlgMigration[2] = { 0.11, 0.83 };
float Progen::cytokineSynthesis[3] = { 1, 2.58, 0 };
float Progen::AggrecanSynth[1] = { 1 };
float Progen::proliferation[1] = {24};
float Progen::differentiation[3] = {0.7, 0.3, 48};

int NP::numOfNP = 0;
float NP::migrationSpeed = 1;    // patch/tick. Not an input parameter
float NP::OCR = 15.31 / 2; // fmol/h/cell divided by 2 for /tick
float NP::apoptosisChance = 5;
float NP::collagenSynthRate = 1;
float NP::aggrecanSynthRate = 1.5;

float NP::CaAlgMigration[2] = { 0.11, 1.30 };
float NP::CollagenSynth[3] = { 10, 6.45, 3.6 };
float NP::AggrecanSynth[3] = { 20, 38, 16.6 };

//DEFAULT CONSTRUCTORS
Cell::Cell() {
	cout << "default cell alloc" << endl;

	// added for debugging: print out what cell type
	//if (typeid(*this) == typeid(stem)) {
	//	cout << "cell type stem";
	//}
	//else if (typeid(*this) == typeid(progen)) {
	//	cout << "cell type progen";
	//}
	//else if (typeid(*this) == typeid(np)) {
	//	cout << "cell type np";
	//}
}

Stem::Stem() {
	cout << "default stem alloc" << endl;
}

Progen::Progen() {
	cout << "default pre-np alloc" << endl;
}

NP::NP() {
	cout << "default np alloc" << endl;
}


//CONSTRUCTORS FOR PATCH POINTERS
Cell::Cell(Patch* patchPtr) {
	this->ix[write_t] = patchPtr->indice[0];
	this->iy[write_t] = patchPtr->indice[1];
	this->iz[write_t] = patchPtr->indice[2];
	this->index[write_t] = patchPtr->index;
	this->alive[write_t] = true;
	this->realDeath[write_t] = false;
	this->doublings[write_t] = 0;


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
	this->life[write_t] = 0;
	
	this->activate[read_t] = false;
	//this->color[read_t] = ccell;
	this->size[read_t] = 2;
	//this->type[read_t] = cell;
}

Stem::Stem(Patch* patchPtr) : Cell(patchPtr) {
	this->alive[write_t] = true;
	this->realDeath[write_t] = false;
	this->doublings[write_t] = 0;
	this->color[write_t] = cstem;
	this->type[write_t] = stem;

	this->doublings[read_t] = 0;
	this->color[read_t] = cstem;
	this->type[read_t] = stem;
}

Progen::Progen(Patch* patchPtr) : Cell(patchPtr) {
	this->alive[write_t] = true;
	this->realDeath[write_t] = false;
	this->doublings[write_t] = 0;
	this->color[write_t] = cprogen;
	this->type[write_t] = progen;

	this->doublings[read_t] = 0;
	this->color[read_t] = cprogen;
	this->type[read_t] = progen;
}

NP::NP(Patch* patchPtr) : Cell(patchPtr) {
	this->alive[write_t] = true;
	this->realDeath[write_t] = false;
	this->doublings[write_t] = 0;
	this->color[write_t] = cnp;
	this->type[write_t] = np;

	this->doublings[read_t] = 0;
	this->color[read_t] = cnp;
	this->type[read_t] = np;

	#ifndef MODEL_SCAFFOLD
		// Unactivated chondrocytes live for 5 to 11 days. 0 corresponds to hours:
		if (Agent::agentWorldPtr->clock == 0) this->life[write_t] = BMWorld::reportTick(0, rand() % 12);
		else this->life[write_t] = BMWorld::reportTick(0, 5 + rand() % 7);
	#else
		//this->life[write_t] = BMWorld::reportTick(0, 5 + rand() % 7);
	this->life[write_t] = 0;
	#endif
}

//CONSTRUCTORS WITH COORDINATES
Cell::Cell(int x, int y, int z) {
	this->ix[write_t] = x;
	this->iy[write_t] = y;
	this->iz[write_t] = z;
	this->index[write_t] = x + y*nx + z*nx*ny;
	this->alive[write_t] = true;
	this->realDeath[write_t] = false;
	this->doublings[write_t] = 0;

	//#ifndef MODEL_SCAFFOLD
	//	// Unactivated chondrocytes live for 5 to 11 days. 0 corresponds to hours.
	//	if (Agent::agentWorldPtr->clock == 0) this->life[write_t] = BMWorld::reportTick(0, rand()%12);
	//	else this->life[write_t] = BMWorld::reportTick(0, 5 + rand()%7);
	//#else
	//	this->life[write_t] = BMWorld::reportTick(0, 5 + rand()%7);
	//#endif

	this->life[write_t] = 0;

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
		this->doublings[write_t] = 0;
	#else
		this->ix[read_t] = x;
		this->iy[read_t] = y;
		this->iz[read_t] = z;
		this->index[read_t] = x + y*nx + z*nx*ny;
		this->alive[read_t] = true;
		this->life[read_t] = 0;
		this->activate[read_t] = false;
		this->color[read_t] = ccell;
		this->size[read_t] = 2;
		this->type[read_t] = cell;
		this->doublings[write_t] = 0;
	#endif
	///* Added by MM to check types of cell stages and add to respective counters: */
	//if (typeid(*this) == typeid(Stem)) {
	//	Stem::numOfStem++;
	//}
	//else if (typeid(*this) == typeid(Progen)) {
	//	Progen::numOfProgen++;
	//}
	//else if (typeid(*this) == typeid(NP)) {
	//	NP::numOfNP++;
	//}
	//Cell::numOfCells++;  
}

Stem::Stem(int x, int y, int z) : Cell(x, y, z) {}

Progen::Progen(int x, int y, int z) : Cell(x, y, z) {}

NP::NP(int x, int y, int z) : Cell(x, y, z) {}

//DESTRUCTORS
Cell::~Cell() {}

Stem::~Stem() {}

Progen::~Progen() {}

NP::~NP() {}

//CELL FUNCTIONS
void Cell::cellFunction() {
	// Calls the individual cell stage functions
	int in = this->index[read_t];

	//Measure mean & patch oxygen 
	float meanO2 = this->meanNeighborConcentration(o2);
	float patchO2 = this->patchChemConcentration(o2, in);

	if (this->alive[read_t] == false) return;
	if (this->alive[read_t] == true && (meanO2 + patchO2) > this->get_OCR()) {
		this->proliferate();
		this->differentiate();
		this->cellSniff();
		this->ecm_synthesis();
		this->cytokine_synthesis();
		this->apoptose();

		// Finally, subtract 'consumed' O2 from current and neighbor patches

		// Location of agent in x,y,z dimensions of world.
		int x = this->ix[read_t];
		int y = this->iy[read_t];
		int z = this->iz[read_t];

		// Number of patches in x,y,z dimensions of world
		int nx = Agent::nx;
		int ny = Agent::ny;
		int nz = Agent::nz;

		float oxDecrease = this->get_OCR()/ 27; // divide OCR across 27 patches (current + neighbors)
		this->addPatchChemSecretion(o2, in, -1 * oxDecrease);
		// Count number of patches of neighbors inside world dimensions:
		for (int dZ = -1; dZ <= 1; dZ++) {
			for (int dY = -1; dY <= 1; dY++) {
				for (int dX = -1; dX <= 1; dX++) {
					if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny || z + dZ < 0 || z + dZ >= nz) continue;
					int in = (x + dX) + (y + dY) * nx + (z + dZ) * nx * ny;
					if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) (this->addPatchChemSecretion(o2, in, -1 * oxDecrease));
				}
			}
		}

		// last thing to do: increase age + 1 tick 
		if (this->life[read_t] >= 0) {
			this->life[write_t] = this->life[read_t] + 1;
		}
	}
}

void Cell::cellSniff() {
	int in = this->index[read_t];

	// gets migration speed from correct hook function
	float speed = get_migration_speed();

	if ((Agent::agentPatchPtr[in]).inDamzone == true) {
		if (rollDice(80) && this->moveTowardChemotaxis() == true){
			#ifdef MODEL_SCAFFOLD
			if (speed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
				// Move up to "migrationSpeed" patches per tick:
				for (int dx = 0; dx < speed; dx++) this->wiggle();
			}
		// If cell is not actively migrating, consider chance of moving to next patch:
		} else if(rollDice(0.25)){		
			this->wiggle(); 
		}
		#else
			this->wiggle();
		#endif

	} else {
		// not in damage zone
#ifdef MODEL_SCAFFOLD
		if (speed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
			for (int dx = 0; dx < speed; dx++) this->wiggle();
		}
		else if (rollDice(0.25)) {
			this->wiggle();
		}
#else
		this->wiggle();
#endif
	}

	// TGF can excite NP cell and overcome gradient: //NOTE MM: double check this for MSCs
	if (can_tgf_excite()) {
#ifdef MODEL_SCAFFOLD
		if (speed > 1 && Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
			for (int dx = 0; dx < speed; dx++) this->wiggle(); // Move up to "migrationSpeed" patches per tick:
		}
		else if (rollDice(0.25)) {
			this->wiggle(); // If cell is not actively migrating, consider chance of moving to next patch: 	
		}
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

void Cell::apoptose() {
#ifdef CALIBRATION
	if (rollDice(get_apoptosis_chance())) {
		this->realDeath[write_t] = true;
		this->die();
		return;
	}
#else
	if (rollDice(1)) { // from Netlogo model
		this->die();
		return;
	}
#endif
}

void Cell::copyAndInitialize(Agent* original, int dx, int dy, int dz) {
	int in = this->index[read_t];

	// Initializes location of new Cell relative to original agent:
	this->ix[write_t] = original->getX() + dx;
	this->iy[write_t] = original->getY() + dy;
	this->iz[write_t] = original->getZ() + dz;
	this->index[write_t] = this->ix[write_t] + this->iy[write_t]*Agent::nx + this->iz[write_t]*Agent::nx*Agent::ny;
  	
	// Initializes new Cell:
	this->alive[read_t] = true;
	this->life[read_t] = 0; // 0 corresponds to ticks
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

void Cell::proliferate() {
	int in = this->index[read_t];
	if (!(Agent::agentPatchPtr[in].type[read_t] == CaAlg)) return; // check for being on a biomaterial patch
	if (!(this->life[read_t] > 0 && this->life[read_t] % Cell::proliferation[0] == 0)) return; // check for 24-hour mark (of the cell's life) to try division
	if (!isProliferative()) return; // check if cell is proliferative; i.e., under the max # of divisions for its type

	// calculating local cytokines
	float meanTNF = this->meanNeighborConcentration(TNF);
	float meanTGF = this->meanNeighborConcentration(TGF);
	float meanIL1 = this->meanNeighborConcentration(IL1beta);

	float prob = get_prolif_prob(meanTGF, meanIL1, meanTNF); // get the proliferation probability for the cell type

	if (rollDice(prob)) {
		this->hatchnewcell(1, this->type[read_t]);
		this->doublings[write_t] = this->doublings[read_t] + 1;
		return;
	}
}

void Cell::differentiate() {
	int in = this->index[read_t];
	if (!(Agent::agentPatchPtr[in].type[read_t] == CaAlg)) return; // check for being on a biomaterial patch
	if (!(this->life[read_t] > 0 && this->life[read_t] % Stem::differentiation[3] == 0)) return; // check for 48-hour mark (of the cell's life) to try differentiation
	if (!isProliferative()) return;

	// calculating local cytokines
	float meanTNF = this->meanNeighborConcentration(TNF);
	float meanTGF = this->meanNeighborConcentration(TGF);
	float meanIL1 = this->meanNeighborConcentration(IL1beta);

	float prob = get_diff_prob(meanTGF, meanIL1, meanTNF);
	int daughterType = get_daughter_type();
	if (daughterType == -1) return; // base Cell has no daughter type; skip 

	if (rollDice(prob)) {
		if (rollDice(Stem::differentiation[0]*100)) { // check for asymmetric differentiation; more likely
			this->hatchnewcell(1, daughterType);
		}
		else { // check for symmetric differentiation; less likely
			Agent::agentPatchPtr[in].clearOccupied();
			Agent::agentPatchPtr[in].occupiedby[write_t] = nothing;

			this->hatchnewcell(1, daughterType, 1); // 'change' cell here to next cell stage
			this->doublings[write_t] = this->doublings[read_t] + 1;
			this->hatchnewcell(1, daughterType); // create new cell in next stage nearby
			this->die(); // 'kill' current cell
		}
	}
}

void Cell::ecm_synthesis() {
	int in = this->index[read_t];

	// Calculates chemical gradients and patch chemical concentrations:
	float meanTNF = this->meanNeighborConcentration(TNF);
	float meanTGF = this->meanNeighborConcentration(TGF);
	float meanIL1 = this->meanNeighborConcentration(IL1beta);
	float patchTNF = this->patchChemConcentration(TNF, in);
	float patchTGF = this->patchChemConcentration(TGF, in);
	float patchIL1beta = this->patchChemConcentration(IL1beta, in);

	// Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];

	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	int neighborCount = 0;
	// Count number of patches of neighbors inside world dimensions:
	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
				if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny || z + dZ < 0 || z + dZ >= nz) continue;
				int in = (x + dX) + (y + dY) * nx + (z + dZ) * nx * ny;
				if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) neighborCount++;
			}
		}
	}

	// Calculate total volume of surrounding patches to check for cytokine thresholds
	float patchVolume = BMWorld::totalVolumeML / (nx * ny * nz);
	//int neighbors = BMWorld::countNeighborPatchType(x, y, z, CaAlg);
	float patchesVolume = patchVolume * neighborCount;

	calculate_ecm_synth_rates(meanTGF, meanIL1, meanTNF, patchesVolume);

	if (fmod(((Agent::agentWorldPtr)->reportHour()), 12) == 0) {
		create_ecm(meanTGF, meanIL1, meanTNF);
	}
}

void Cell::cytokine_synthesis() {
	int in = this->index[read_t];

	// Calculates patch chemical concentrations:
	float patchTNF = this->patchChemConcentration(TNF, in);
	float patchTGF = this->patchChemConcentration(TGF, in);
	float patchIL1beta = this->patchChemConcentration(IL1beta, in);

	create_cytokines(patchTGF, patchIL1beta, patchTNF);
}

void Cell::calculate_ecm_synth_rates(float meanTGF, float meanIL1, float meanTNF, float patchesVolume) {}

void Cell::create_ecm(float meanTGF, float meanIL1, float meanTNF) {}

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

		// Move to new patch and sprout oaggrecan
		in = (x + dx) + (y + dy) * nx + (z + dz) * nx * ny;
		this->move(dx, dy, dz, read_index);

		Agent::agentECMPtr[in].oaggrecan[write_t] = Agent::agentECMPtr[in].oaggrecan[read_t] + 1 + rand() % 2;
#ifdef OPT_ECM
		Agent::agentECMPtr[in].set_dirty();
#endif
	}
}

void Cell::create_cytokines(float patchTGF, float patchIL1beta, float patchTNF) {}

bool Cell::isProliferative() {
	return this->doublings[read_t] < get_max_doublings();
}

int Cell::get_max_doublings() { return 50; } //base default

float Cell::get_prolif_prob(float meanTGF,
	float meanIL1,
	float meanTNF) { return 10; } //base default

int Cell::get_daughter_type() { return -1; } // sentinel 'no type'

float Cell::get_diff_prob(float meanTGF,
	float meanIL1,
	float meanTNF) { return 5; } // base default

float Cell::get_migration_speed() { return 0; }

float Cell::get_apoptosis_chance() { return 1; }

float Cell::get_OCR() { return 0; }

void Cell::hatchnewcell(int number, int agentType, int here) {
	int newcells = 0;
	int lx = 0;
	int ly = 0;
	int lz = 0;
	int in = 0;

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
		if (here == 0) { // Default option; hatching on number of neighboring patches
			// Distance away from target neighboring patch in x,y,z dimensions
			int dx = Agent::dX[neighbor[i]];
			int dy = Agent::dY[neighbor[i]];
			int dz = Agent::dZ[neighbor[i]];

			// Patch row major index of target neighboring patch:
			in = (x + dx) + (y + dy) * nx + (z + dz) * nx * ny;

			// Hatching coordinates:
			int lx = x + dx;
			int ly = y + dy;
			int lz = z + dz;

			// Try a new target neighboring patch if this one is not inside the world dimensions, or is occupied.
			if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny || z + dz < 0 || z + dz >= nz) continue;
			int targetType = agentPatchPtr[in].type[read_t];
		}
		else { // here == 1; option to hatch a new cell on current patch
			in = this->getIndex();

			// Hatching coordinates:
			int lx = x;
			int ly = y;
			int lz = z;
		}
		// Create a new cell of agentType at the valid target neighboring patch:
		Cell* newcell = nullptr;
		switch (agentType) {
		case stem:
		{
			newcell = new Stem(lx, ly, lz);
		}
		break;
		case progen:
		{
			newcell = new Progen(lx, ly, lz);
		}
		break;
		case np:
		{
			newcell = new NP(lx, ly, lz);
		}
		break;
		}
		newcells++;

		// Update target neighboring patch as occupied:
		Agent::agentPatchPtr[in].setOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = agentType;

		/* If executing OMP version, add the pointer to this new cell to the thread-local list first.
	* BMWorld::UpdateCells() will take care of putting it in the global list at the end
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

/* -------------------------------------------------------------------------- */
/*                                    STEM                                    */
/* -------------------------------------------------------------------------- */

int Stem::get_max_doublings() { return 100; }

float Stem::get_prolif_prob(float meanTGF,
	float meanIL1,
	float meanTNF) {

	int TGFrelated = 0;

#ifdef CALIBRATION
	if (meanTGF <= Stem::proliferation[0]) {
#else  
	if (meanTGF <= 10) {
#endif
		TGFrelated = 1;  // Low TGF (0.1-1ng) stimulate proliferation and attraction. 
	}
	else {
		TGFrelated = -1; // High TGF (1-10ng) inhibits proliferation. 
	}

#ifdef CALIBRATION
	//float prolif = log10(1 - Stem::proliferation[1] * meanTNF - Stem::proliferation[2] * meanIL1 + TGFrelated * meanTGF - Stem::proliferation[3] * BMWorld::E);
	float prolif = 40; //testing
#else  
	float prolif = log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF);
#endif  
	
	return prolif;
}

int Stem::get_daughter_type() { return progen; }

float Stem::get_diff_prob(float meanTGF,
	float meanIL1,
	float meanTNF) {

	//return (Stem::differentiation[1]*100) + (Stem::differentiation[2] * meanTGF);
	return 20;
}

void Stem::calculate_ecm_synth_rates(float meanTGF, float meanIL1, float meanTNF, float patchesVolume) {
#ifdef CALIBRATION
	Stem::collagenSynthRate = Stem::CollagenSynth[0]
		+ (log10(1 + meanTGF) / (1 + meanTNF + meanIL1));

	if (meanTGF < (Stem::AggrecanSynth[0] / patchesVolume)) {
		Stem::aggrecanSynthRate = Stem::collagenSynthRate / 1.2;
	}
	else {
		Stem::aggrecanSynthRate = Stem::collagenSynthRate * 1.2;
	}
#else
	Stem::collagenSynthRate = 10
		+ (log10(1 + meanTGF) / (1 + meanTNF + meanIL1));

	if (meanTGF < 100000) {
		Stem::aggrecanSynthRate = Stem::collagenSynthRate / 1.2;
	}
	else {
		Stem::aggrecanSynthRate = Stem::collagenSynthRate * 1.2;
	}
#endif
}

void Stem::create_ecm(float meanTGF, float meanIL1, float meanTNF) {
	int in = this->index[read_t];

#ifdef MODEL_SCAFFOLD
	if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
		// loops based on synth rates calculated in calculate_ecm_synth_rates()
		for (int i = 0; i < Stem::collagenSynthRate; i++)
			this->makeOCollagen(meanTGF, meanIL1);
		for (int i = 0; i < Stem::aggrecanSynthRate; i++)
			this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
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

void Stem::create_cytokines(float patchTGF, float patchIL1beta, float patchTNF) {
	int in = this->index[read_t];
	// Change in chemicals due to cells:
#ifdef CALIBRATION
	this->addPatchChemSecretion(TGF, in, Stem::cytokineSynthesis[0] + Cell::cytokineSynthesis[1] * (patchTGF)+Cell::cytokineSynthesis[2] * (patchIL1beta)+Cell::cytokineSynthesis[3] * (patchTNF));//this->addPatchChemSecretion(TGF, in,  Chondrocyte::cytokineSynthesis[0] + Chondrocyte::cytokineSynthesis[1]*(1 + Chondrocyte::cytokineSynthesis[2]*patchTNF);
	//this->addPatchChemSecretion(TGF, in, 0; //DEBUG : constant TGF
	this->addPatchChemSecretion(TNF, in, Stem::cytokineSynthesis[1] + (Cell::cytokineSynthesis[5] * ((patchIL1beta) / (1 + Cell::cytokineSynthesis[6] * patchTGF))));//this->addPatchChemSecretion(TNF, in, Chondrocyte::cytokineSynthesis[3] + Chondrocyte::cytokineSynthesis[4]/(1 + patchTGF*Chondrocyte::cytokineSynthesis[5]);
	//this->addPatchChemSecretion(TNF, in, 0; //DEBUG : constant TNF
	this->addPatchChemSecretion(IL1beta, in, Stem::cytokineSynthesis[2] + (Cell::cytokineSynthesis[8] * ((patchTNF) / (1 + Cell::cytokineSynthesis[9] * patchTGF))));//this->addPatchChemSecretion(IL1beta, in, Chondrocyte::cytokineSynthesis[6] + (Chondrocyte::cytokineSynthesis[7]*patchTNF)/(Chondrocyte::cytokineSynthesis[8] + Chondrocyte::cytokineSynthesis[9]*patchTGF);
#else
	this->addPatchChemSecretion(TGF, in, 5 + (2.25 * patchTGF + 1.3 * patchIL1beta + 5.11 * patchTNF));//9.98 + 2.58*patchTGF + 5.11*patchTNF;				//2.11 + 3.7*patchTGF;
	this->addPatchChemSecretion(TNF, in, 0 + (2.42 * patchIL1beta) / (1 + 4.22 * patchTGF));//5.16 + (2.42*patchIL1beta)/(1 + 4.22*patchTGF);	//2.4*patchIL1beta + 4.8/(1 + 1.27*patchTGF);		
	this->addPatchChemSecretion(IL1beta, in, 0 + (5.43 * patchTNF) / (1 + 3.26 * patchTGF));//2.11 + (5.43*patchTNF)/(1 + 3.26*patchTGF);		//4;
#endif
}

float Stem::get_migration_speed() {
	if (BMWorld::clock == 0) {
#ifdef CALIBRATION
		float migration_ummin = Stem::CaAlgMigration[0] * log(BMWorld::E) + Stem::CaAlgMigration[1]; // um/min

		if (rollDice(0.5)) {  // Convert migration speed in um/min to patches/tick where default patchlength is 10um and default tick is 30 min
			Stem::migrationSpeed = ceil(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
		}
		else {
			Stem::migrationSpeed = floor(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
		}
#else
		float migration_ummin = 0.1096 * log(BMWorld::E) + 0.35; // um/min //float migration_ummin =  0.1213*log10(Agent::agentWorldPtr->E) + 0.223; // um/min

		if (rollDice(0.5)) {  // Convert migration speed in um/min to patches/tick where default patchlength is 10um and default tick is 30 min
			Stem::migrationSpeed = ceil(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
		}
		else {
			Stem::migrationSpeed = floor(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
		}
#endif
		cout << "        Stem cell migration Speed (patch/tick) = " << Stem::migrationSpeed << endl;
	}
	return Stem::migrationSpeed;
}

float Stem::get_apoptosis_chance() {
	return Stem::apoptosisChance;
}

float Stem::get_OCR() {
	return Stem::OCR;
}

/* -------------------------------------------------------------------------- */
/*                                PROGENITOR                                  */
/* -------------------------------------------------------------------------- */
int Progen::get_max_doublings() { return 65; }

float Progen::get_prolif_prob(float meanTGF,
	float meanIL1,
	float meanTNF) {

#ifdef CALIBRATION
	//float prolif = log10(1 + meanTNF - meanIL1 + meanTGF);
	float prolif = 40; //testing
#else  
	float prolif = log10(1 + meanTNF - meanIL1 + meanTGF);
#endif 

	return prolif;
}

int Progen::get_daughter_type() { return np; }

float Progen::get_diff_prob(float meanTGF,
	float meanIL1,
	float meanTNF) {

	return 20;
}

void Progen::calculate_ecm_synth_rates(float meanTGF, float meanIL1, float meanTNF, float patchesVolume) {
#ifdef CALIBRATION
	Progen::aggrecanSynthRate = Progen::AggrecanSynth[0]
		+ (log10(1 + meanTGF) / (1 + meanTNF + meanIL1));
#else
	Progen::aggrecanSynthRate = Progen::AggrecanSynth[0];
#endif
}

void Progen::create_ecm(float meanTGF, float meanIL1, float meanTNF) {
	int in = this->index[read_t];

#ifdef MODEL_SCAFFOLD
	if (Agent::agentPatchPtr[in].type[read_t] == CaAlg) {
		// progen on scaffold only loops aggrecan, never collagen
		for (int i = 0; i < Progen::aggrecanSynthRate; i++)
			this->makeOAggrecan(meanTNF, meanTGF, meanIL1);
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

void Progen::create_cytokines(float patchTGF, float patchIL1beta, float patchTNF) {
	int in = this->index[read_t];
	// Change in chemicals due to cells:
#ifdef CALIBRATION
	this->addPatchChemSecretion(TGF, in, Progen::cytokineSynthesis[0] + Cell::cytokineSynthesis[1] * (patchTGF)+Cell::cytokineSynthesis[2] * (patchIL1beta)+Cell::cytokineSynthesis[3] * (patchTNF));//this->addPatchChemSecretion(TGF, in,  Chondrocyte::cytokineSynthesis[0] + Chondrocyte::cytokineSynthesis[1]*(1 + Chondrocyte::cytokineSynthesis[2]*patchTNF);
	//this->addPatchChemSecretion(TGF, in, 0; //DEBUG : constant TGF
	this->addPatchChemSecretion(TNF, in, Progen::cytokineSynthesis[1] + (Cell::cytokineSynthesis[5] * ((patchIL1beta) / (1 + Cell::cytokineSynthesis[6] * patchTGF))));//this->addPatchChemSecretion(TNF, in, Chondrocyte::cytokineSynthesis[3] + Chondrocyte::cytokineSynthesis[4]/(1 + patchTGF*Chondrocyte::cytokineSynthesis[5]);
	//this->addPatchChemSecretion(TNF, in, 0; //DEBUG : constant TNF
	this->addPatchChemSecretion(IL1beta, in, Progen::cytokineSynthesis[2] + (Cell::cytokineSynthesis[8] * ((patchTNF) / (1 + Cell::cytokineSynthesis[9] * patchTGF))));//this->addPatchChemSecretion(IL1beta, in, Chondrocyte::cytokineSynthesis[6] + (Chondrocyte::cytokineSynthesis[7]*patchTNF)/(Chondrocyte::cytokineSynthesis[8] + Chondrocyte::cytokineSynthesis[9]*patchTGF);
#else
	this->addPatchChemSecretion(TGF, in, 1 + (0.25 * patchTGF + 1.3 * patchIL1beta + 5.11 * patchTNF));//9.98 + 2.58*patchTGF + 5.11*patchTNF;				//2.11 + 3.7*patchTGF;
	this->addPatchChemSecretion(TNF, in, 2.58 + (2.42 * patchIL1beta) / (1 + 4.22 * patchTGF));//5.16 + (2.42*patchIL1beta)/(1 + 4.22*patchTGF);	//2.4*patchIL1beta + 4.8/(1 + 1.27*patchTGF);		
	this->addPatchChemSecretion(IL1beta, in, 0 + (5.43 * patchTNF) / (1 + 3.26 * patchTGF));//2.11 + (5.43*patchTNF)/(1 + 3.26*patchTGF);		//4;
#endif
}

float Progen::get_migration_speed() {
	if (BMWorld::clock == 0) {
#ifdef CALIBRATION
		float migration_ummin = Progen::CaAlgMigration[0] * log(BMWorld::E) + Progen::CaAlgMigration[1]; // um/min

		if (rollDice(0.5)) {  // Convert migration speed in um/min to patches/tick where default patchlength is 10um and default tick is 30 min
			Progen::migrationSpeed = ceil(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
		}
		else {
			Progen::migrationSpeed = floor(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
		}
#else
		float migration_ummin = 0.1096 * log(BMWorld::E) + ((NP::migrationSpeed - Stem::migrationSpeed) / 2); // um/min //float migration_ummin =  0.1213*log10(Agent::agentWorldPtr->E) + 0.223; // um/min

		if (rollDice(0.5)) {  // Convert migration speed in um/min to patches/tick where default patchlength is 10um and default tick is 30 min
			Progen::migrationSpeed = ceil(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
		}
		else {
			Progen::migrationSpeed = floor(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
		}
#endif
		cout << "        Progenitor cell migration Speed (patch/tick) = " << Progen::migrationSpeed << endl;
	}
	return Progen::migrationSpeed;
}

float Progen::get_apoptosis_chance() {
	return Progen::apoptosisChance;
}

float Progen::get_OCR() {
	return Progen::OCR;
}

/* -------------------------------------------------------------------------- */
/*                                     NP                                     */
/* -------------------------------------------------------------------------- */
int NP::get_max_doublings() { return 27; }

float NP::get_prolif_prob(float meanTGF,
	float meanIL1,
	float meanTNF) {

	int TGFrelated = 0;

#ifdef CALIBRATION
	if (meanTGF <= Cell::proliferation[1]) {
#else
	if (meanTGF <= 10) {
#endif
		TGFrelated = 1;  // Low TGF (0.1-1 ng) stimulate chond proliferation and attraction.
	}
	else {
		TGFrelated = -1;  // High TGF (1-10 ng) inhibits proliferation.
	}

#ifdef CALIBRATION
	float prolif = Cell::proliferation[2] * (log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF)) + Cell::proliferation[3];
#else  
	float prolif = log10(1 + meanTNF + meanIL1 + TGFrelated * meanTGF);
#endif 

	return prolif;
}

void NP::calculate_ecm_synth_rates(float meanTGF, float meanIL1, float meanTNF, float patchesVolume) {
#ifdef CALIBRATION
	NP::collagenSynthRate = NP::CollagenSynth[0]
		* (NP::CollagenSynth[1] * BMWorld::reportDay() + NP::CollagenSynth[2]);
	NP::aggrecanSynthRate = NP::AggrecanSynth[0]
		* (NP::AggrecanSynth[1] * BMWorld::reportDay() + NP::AggrecanSynth[2]);
#else
	NP::collagenSynthRate = 10 * (6.45 * BMWorld::reportDay() + 3.6);
	NP::aggrecanSynthRate = 20 * (38 * BMWorld::reportDay() + 16.6);
#endif
}

void NP::create_ecm(float meanTGF, float meanIL1, float meanTNF) {
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
}

void NP::create_cytokines(float patchTGF, float patchIL1beta, float patchTNF) {
	int in = this->index[read_t];
	// Change in chemicals due to cells:
#ifdef CALIBRATION
	this->addPatchChemSecretion(TGF, in, Cell::cytokineSynthesis[0] + Cell::cytokineSynthesis[1] * (patchTGF)+Cell::cytokineSynthesis[2] * (patchIL1beta)+Cell::cytokineSynthesis[3] * (patchTNF));//this->addPatchChemSecretion(TGF, in,  Chondrocyte::cytokineSynthesis[0] + Chondrocyte::cytokineSynthesis[1]*(1 + Chondrocyte::cytokineSynthesis[2]*patchTNF);
	//this->addPatchChemSecretion(TGF, in, 0; //DEBUG : constant TGF
	this->addPatchChemSecretion(TNF, in, Cell::cytokineSynthesis[4] + (Cell::cytokineSynthesis[5] * ((patchIL1beta) / (1 + Cell::cytokineSynthesis[6] * patchTGF))));//this->addPatchChemSecretion(TNF, in, Chondrocyte::cytokineSynthesis[3] + Chondrocyte::cytokineSynthesis[4]/(1 + patchTGF*Chondrocyte::cytokineSynthesis[5]);
	//this->addPatchChemSecretion(TNF, in, 0; //DEBUG : constant TNF
	this->addPatchChemSecretion(IL1beta, in, Cell::cytokineSynthesis[7] + (Cell::cytokineSynthesis[8] * ((patchTNF) / (1 + Cell::cytokineSynthesis[9] * patchTGF))));//this->addPatchChemSecretion(IL1beta, in, Chondrocyte::cytokineSynthesis[6] + (Chondrocyte::cytokineSynthesis[7]*patchTNF)/(Chondrocyte::cytokineSynthesis[8] + Chondrocyte::cytokineSynthesis[9]*patchTGF);
#else
	this->addPatchChemSecretion(TGF, in, 10 + 0.05 * (patchTGF + 10 * patchTNF));//9.98 + 2.58*patchTGF + 5.11*patchTNF;				//2.11 + 3.7*patchTGF;
	this->addPatchChemSecretion(TNF, in, 5 + (2.4 * patchIL1beta) / (1 + 4 * patchTGF));//5.16 + (2.42*patchIL1beta)/(1 + 4.22*patchTGF);	//2.4*patchIL1beta + 4.8/(1 + 1.27*patchTGF);		
	this->addPatchChemSecretion(IL1beta, in, 2 + (5 * patchTNF) / (1 + 3.2 * patchTGF));//2.11 + (5.43*patchTNF)/(1 + 3.26*patchTGF);		//4;
#endif
}

float NP::get_migration_speed() {
#ifdef CALIBRATION
	float migration_ummin = 0.1096 * log(BMWorld::E) + 0.2431; // um/min

	if (rollDice(0.5)) {  // Convert migration speed in um/min to patches/tick where default patchlength is 10um and default tick is 30 min
		NP::migrationSpeed = ceil(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
	}
	else {
		NP::migrationSpeed = floor(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
	}
#else
	float migration_ummin = 0.1096 * log(BMWorld::E) + 0.2431; // um/min //float migration_ummin =  0.1213*log10(Agent::agentWorldPtr->E) + 0.223; // um/min

	if (rollDice(0.5)) {  // Convert migration speed in um/min to patches/tick where default patchlength is 10um and default tick is 30 min
		NP::migrationSpeed = ceil(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
	}
	else {
		NP::migrationSpeed = floor(migration_ummin * 30 / (Agent::agentWorldPtr->patchlength * pow(10, 3)));    //patch/tick 
	}
#endif
	cout << "        NP cell migration Speed (patch/tick) = " << NP::migrationSpeed << endl;
	return NP::migrationSpeed;
}

bool NP::can_tgf_excite() {
	return this->meanNeighborConcentration(TGF) > 0;
}

float NP::get_apoptosis_chance() {
	return NP::apoptosisChance;
}

float NP::get_OCR() {
	return NP::OCR;
}