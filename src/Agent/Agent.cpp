/* 
 * File: Agent.cpp
 * 
 * File Contents: Contains the Agent class.
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 *               Meghana Munipalle
 */

#include "Agent.h"
#include "../Chemistry/chemical_environment.h"
#include "../World/Usr_World/biomaterialWorld.h"
#include "../enums.h"
#include <iostream>
#include <vector>
#include <tgmath.h>
#include <iomanip> // Include for setprecision and fixed


BMWorld* Agent::agentWorldPtr = NULL;
Patch* Agent::agentPatchPtr = NULL;
ECM* Agent::agentECMPtr = NULL;

ChemicalEnvironment* Agent::chemicalEnvironment()
{
	if (Agent::agentWorldPtr)
		return Agent::agentWorldPtr->chemical_environment();
	return nullptr;
}

int Agent::nx = 0;
int Agent::ny = 0;
int Agent::nz = 0;
int Agent::dX[27] = {-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1,  0, 1, -1, 0, 1, -1, 0, 1};
int Agent::dY[27] = {-1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1};
int Agent::dZ[27] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};

#ifdef MODEL_3D
	int Agent::neighbor[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
#else
	int Agent::neighbor[8] = {9, 10, 11, 12, 14, 15, 16, 17}; // We do not include neighbor 13 (no movement to self (0,0,0))
#endif

using namespace std;

Agent::Agent() {}
Agent::~Agent() {}

void Agent::openRngStream(int birth_index) {
	this->rng = abm::rng::Stream(static_cast<uint32_t>(birth_index),
	                             static_cast<uint32_t>(World::clock));
}

bool Agent::isAlive() {
	if (this->life[write_t] < 0) this->alive[write_t] = false;
	else this->alive[write_t] = true;
	return this->alive[read_t];
} 

bool Agent::isRealDead() {
	return this->realDeath[read_t];
}

void Agent::cellFunction() {}
void Agent::copyAndInitialize(Agent* original, int dx, int dy, int dz) {}
void die() {}

int Agent::getX() {
	return this->ix[read_t];
}

int Agent::getY() {
	return this->iy[read_t];
}

int Agent::getZ() {
	return this->iz[read_t];
}

int Agent::getIndex() {
	return this->index[read_t];
}

bool Agent::rollDice(float percent) {
	return abm::rng::roll_percent(this->rng, percent);
}

bool Agent::move(int dX, int dY, int dZ, int read_index) {
  // Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];

  // Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

  // Location of target patch in x,y,z dimensions of world.
	int targetX = x + dX;
	int targetY = y + dY;
	int targetZ = z + dZ;

   /* Abort movement if trying to move to own patch, or to a patch outside world dimensions, 
	* or to the patch the agent was at at the beginning of the tick, or to an occupied patch. */
	if (dX == 0 && dY== 0 && dZ == 0) return false;
	if (targetX < 0 || targetX >= nx || targetY < 0 || targetY >= ny || targetZ < 0 || targetZ >= nz) return false;
	
	int newIndex = targetX + targetY*nx + targetZ*nx*ny;
	if (newIndex == this->index[read_t]) return false;
  
   /* If patch is unoccupied, move current instance to new patch at (x + dX, y + dY, z + dZ). 
	* setOccupied() sets the new patch as occupied and returns true if the patch was already occupied 
    */
	if (!Agent::agentPatchPtr[newIndex].setOccupied()) { 
		int in = this->index[read_index]; // Get the current location
		
		// Update residing patch to 'available'
		Agent::agentPatchPtr[in].clearOccupied();
		Agent::agentPatchPtr[in].occupiedby[write_t] = unoccupied;
		
		// Update location coordinates
		this->ix[write_t] = targetX;
		this->iy[write_t] = targetY;
		this->iz[write_t] = targetZ;
		this->index[write_t] = newIndex;
		Agent::agentPatchPtr[newIndex].occupiedby[write_t] = this->type[read_t];
		return true;
	} else {
		//cout << "ERROR error in cell moving" << endl;
		return false;
	}
}

void Agent::wiggle() {       
	// Check if the location has been modified in this tick:                                               
	int read_index;
	if (isModified(this->index)) read_index = write_t;	// If it has, work off of the intermediate value
	else read_index = read_t;    						// If it has NOT, work off of the original value
	
	// Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_index];
	int y = this->iy[read_index];
	int z = this->iz[read_index];
	int currentindex = this->index[read_index];
	
	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;
	
	int trial = 0;
	do {
		// Pick a neighbor to move to at random:
		#ifndef MODEL_3D
			int i = abm::rng::uniform_int(this->rng, 8);
		#else
			int i = abm::rng::uniform_int(this->rng, 27);
		#endif
			int dx = Agent::dX[neighbor[i]];
			int dy = Agent::dY[neighbor[i]];
			int dz = Agent::dZ[neighbor[i]];

		// Calculate target index:
		int newindex = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

		// If the z-direction movement is invalid, pick a new neighbor
		if (z + dz < 0 || z + dz >= nz) continue;

		// If trying to move off the side boundaries, die: //NOTE by MM: why? why not just count it as an invalid movement?
		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny) {
			//this->life[write_t] = 0;
			//this->die();
			continue;
		}

		// If the target patch is occupied, pick a new neighbor
		if (Agent::agentPatchPtr[newindex].isOccupiedWrite()) continue;
	
		if (Agent::agentPatchPtr[newindex].type[read_t] == CaAlg) {
			dx = -dx;
			dy = -dy;

			vector<int> xtarget;
			vector<int> ytarget;
			vector<int> ztarget;

			// Look for neighbors that are not capillary that are inside world dimensions:
			for (int dxx = -1; dxx <= 1; dxx++) {
				for (int dyy = -1; dyy <= 1; dyy++) {
					for (int dzz = -1; dzz <= 1; dzz++) {
						if (x + dxx < 0 || x + dxx >= nx || y + dyy < 0 || y + dyy >= ny || z + dzz < 0 || z + dzz >= nz) continue;
						
						int in = (x + dxx) + (y + dyy)*nx + (z + dzz)*nx*ny;
						xtarget.push_back(dxx);
						ytarget.push_back(dyy);
						ztarget.push_back(dzz);
					}
				}
			}

			// Move to a random neighbor that is not capillary that is inside world dimensions:
			int randInt = abm::rng::uniform_int(this->rng, xtarget.size());
			dx = xtarget[randInt];
			dy = ytarget[randInt];
			dz = ztarget[randInt];

		} else {
			cout << "exception! encountered by " << this->index[read_t] << " " << Agent::agentPatchPtr[newindex].type[read_t] <<  " " << Agent::agentPatchPtr[this->index[read_t]].type[read_t] << endl;
		}

		if (this->move(dx, dy, dz, read_index) == true) break;  // If move() was successful, get out of the while loop
		else continue;  										// If move() was NOT successful, pick a new neighbor
		
	} while(++trial < 8);

} // End Agent::wiggle()

float Agent::patchChemConcentration(SpeciesId species, int patch_index) {
	const ChemicalEnvironment* env = Agent::chemicalEnvironment();
	if (env)
		return env->concentration_at(patch_index, species);
	return 0.f;
}

void Agent::addPatchChemSecretion(SpeciesId species, int patch_index, float delta) {
	ChemicalEnvironment* env = Agent::chemicalEnvironment();
	if (env)
		env->accumulate_secretion(patch_index, species, delta);
}

float Agent::patchChemotaxis(int patch_index) {
	const ChemicalEnvironment* env = Agent::chemicalEnvironment();
	if (env)
		return env->chemotaxis_at(patch_index);
	return 0.f;
}

float Agent::meanNeighborConcentration(SpeciesId species) {
	float totalchemical = 0.f;
	int numberofpatches = 0;

	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];

	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
				if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny ||
				    z + dZ < 0 || z + dZ >= nz)
					continue;
				int in = (x + dX) + (y + dY) * nx + (z + dZ) * nx * ny;
				totalchemical += this->patchChemConcentration(species, in);
				numberofpatches++;
			}
		}
	}
	return totalchemical / numberofpatches;
}

int Agent::countNeighborECM(int ECMIndex) {
	int numberofecm = 0;
	int numberofColl = 0;
	int numberofAgg = 0;
  	
	// Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];
  	
	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

   // Count number of ECM proteins of type ECMIndex in all neighbors inside world dimensions:
	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
				if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny || z + dZ < 0 || z + dZ >= nz) continue;
				int in = (x + dX) + (y + dY)*nx + (z + dZ)*nx*ny;
				switch (ECMIndex) {
					case orig_coll:
						numberofecm += Agent::agentECMPtr[in].ocollagen[read_t];
						break;
					case new_coll:
						numberofecm += Agent::agentECMPtr[in].ncollagen[read_t];
						break;
					case frag_coll:
						numberofecm += Agent::agentECMPtr[in].fcollagen[read_t];
						break;

					case orig_agg:
						numberofecm += Agent::agentECMPtr[in].oaggrecan[read_t];
						break;
					case new_agg:
						numberofecm += Agent::agentECMPtr[in].naggrecan[read_t];
						break;
					case frag_agg:
						numberofecm += Agent::agentECMPtr[in].faggrecan[read_t];
						break;
				}
			}
		}
	}
	return numberofecm;
}

int Agent::countNeighborCells(int cellIndex) {
	int numberofcells = 0;
  	
	// Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];
  	
	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

    // Count number of cells of type cellIndex in all neighbors inside world dimensions:
	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
				if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny || z + dZ < 0 || z + dZ >= nz) continue;
				
				int in = (x + dX) + (y + dY)*nx + (z + dZ)*nx*ny;
				if (Agent::agentPatchPtr[in].isOccupied() == false) continue;
				if (Agent::agentPatchPtr[in].occupiedby[read_t] == cellIndex) numberofcells++;
			}
		}
	}
	
	return numberofcells;
}

bool Agent::moveTowardChemotaxis() {
	int read_index;

	if (isModified(this->index))
		read_index = write_t;
	else
		read_index = read_t;

	int ix = this->ix[read_index];
	int iy = this->iy[read_index];
	int iz = this->iz[read_index];
	int index = this->index[read_index];

	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	double highestchem = this->patchChemotaxis(index);
	int dx = 0, dy = 0, dz = 0;

	#ifdef MODEL_SCAFFOLD
	int radius;
		switch (this->type[read_t]) {
			case stem: {
				int radius = Stem::migrationSpeed;
			}
			case progen: {
				int radius = Progen::migrationSpeed;
			}
			case np: {
				int radius = NP::migrationSpeed;
			}
		}

		for (int deltaz = -radius; deltaz <= radius; deltaz++) {
			for (int deltay = -radius; deltay <= radius; deltay++) {
				for (int deltax = -radius; deltax <= radius; deltax++) {
					if (ix + deltax < 0 || ix + deltax >= nx || iy + deltay < 0 || iy + deltay >= ny || iz + deltaz < 0 || iz + deltaz >= nz) continue;

					int in = (ix + deltax) + (iy + deltay)*nx + (iz + deltaz)*nx*ny;
					const float neighbor = this->patchChemotaxis(in);
					if (neighbor > highestchem) {
						highestchem = neighbor;
						dx = deltax;
						dy = deltay;
						dz = deltaz;
					}
				}
			}
		}

	#else
		for (int deltaz = -1; deltaz <= 1; deltaz++) {
			for (int deltay = -1; deltay <= 1; deltay++) {
				for (int deltax = -1; deltax <= 1; deltax++) {
					if (ix + deltax < 0 || ix + deltax >= nx || iy + deltay < 0 || iy + deltay >= ny || iz + deltaz < 0 || iz + deltaz >= nz) continue;

					int in = (ix + deltax) + (iy + deltay)*nx + (iz + deltaz)*nx*ny;
					const float neighbor = this->patchChemotaxis(in);
					if (neighbor > highestchem) {
						highestchem = neighbor;
						dx = deltax;
						dy = deltay;
						dz = deltaz;
					}
				}
			}
		}
	#endif

	if (dx == 0 && dy == 0 && dz == 0) return false;
	int newIndex = (ix + dx) + (iy + dy)*nx + (iz + dz)*nx*ny;
	return this->move(dx, dy, dz, read_index);
}

void Agent::updateAgent() {
	this->ix[read_t] = this->ix[write_t];
	this->iy[read_t] = this->iy[write_t];
	this->iz[read_t] = this->iz[write_t];
	this->index[read_t] = this->index[write_t];
	this->alive[read_t] = this->alive[write_t];
	this->life[read_t]  =  this->life[write_t];
	this->activate[read_t] = this->activate[write_t];
	this->color[read_t] = this->color[write_t];
	this->size[read_t]  = this->size[write_t];
	this->type[read_t]  = this->type[write_t];
	this->realDeath[read_t] = this->realDeath[write_t];
	this->doublings[read_t] = this->doublings[write_t];
}