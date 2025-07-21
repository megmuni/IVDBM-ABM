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
#include "../World/Usr_World/woundHealingWorld.h"
#include "../enums.h"
#include <iostream>
#include <vector>
#include <tgmath.h>
#include <iomanip> // Include for setprecision and fixed


WHWorld* Agent::agentWorldPtr = NULL;
Patch* Agent::agentPatchPtr = NULL;
ECM* Agent::agentECMPtr = NULL; 

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

float Agent::viabilityRate = 97; // %
float Agent::proliferationRate = 27;
float Agent::collagenSynthRate = 1; 
float Agent::aggrecanSynthRate = 1.5;//1;
float Agent::HASynthRate = 0;
bool  Agent::CaAlgFlag = false; 

float Agent::CaAlgViability[3] = {0.7321, 97.452, 6};
float Agent::CaAlgProlif[5] = {13.06, 3.69, 101.18, 0.056, 30.44}; //{12, 11.8, 21.5, 20};
float Agent::CollagenSynth[3] = {10, 6.45, 3.6}; //{8922, 255400, 0.02429};
float Agent::AggrecanSynth[3] = {10, 38, 16.6}; //{183.4, 4702, 0.0009759};
float Agent::HASynth[3] = {0,0,0};

using namespace std;

Agent::Agent() {}
Agent::~Agent() {}

bool Agent::isAlive() {
	if (this->life[write_t] <= 0) this->alive[write_t] = false;
	else this->alive[write_t] = true;
	return this->alive[read_t];
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
	int tid = 0;
		#ifdef _OMP
			tid = omp_get_thread_num();	// Get thread id in order to access the seed that belongs to this thread
		#endif
	int randNum = rand_r(&(agentWorldPtr->seeds[tid]))%100;
	if (randNum < percent) return 1;
	else return 0;
}

#ifdef MODEL_SCAFFOLD

	/* ----------------------------- 		MIGRATION SPEED 		 ---------------------------- */
	void Agent::calculateMigrationSpeed(){
		/* Calculate migration speed (patch/tick) of Chondrocyte in Ca-Alg Gel given elastic modulus (E)
		*        
		*        Cells move through gel in brownian motion and along chemical gradients
		*        Cells favor attachment and move faster in hydrogels with high elastic modulus
		*/
		#ifdef CALIBRATION
			float migration_ummin = 0.1096*log(Agent::agentWorldPtr->E) + 0.2431; // um/min

			if (rollDice(0.5)){  // Convert migration speed in um/min to patches/tick where default patchlength is 10um and default tick is 30 min
				Agent::migrationSpeed = ceil(migration_ummin*30/(Agent::agentWorldPtr->patchlength*pow(10,3)));    //patch/tick 
			} else {
				Agent::migrationSpeed = floor(migration_ummin*30/(Agent::agentWorldPtr->patchlength*pow(10,3)));    //patch/tick 
			}
		#else
			float migration_ummin = 0.1096*log(Agent::agentWorldPtr->E) + 0.2431; // um/min //float migration_ummin =  0.1213*log10(Agent::agentWorldPtr->E) + 0.223; // um/min

			if (rollDice(0.5)){  // Convert migration speed in um/min to patches/tick where default patchlength is 10um and default tick is 30 min
				Agent::migrationSpeed = ceil(migration_ummin*30/(Agent::agentWorldPtr->patchlength*pow(10,3)));    //patch/tick 
			} else {
				Agent::migrationSpeed = floor(migration_ummin*30/(Agent::agentWorldPtr->patchlength*pow(10,3)));    //patch/tick 
			}
		#endif
		
		cout << "        Migration Speed (patch/tick) = " << Agent::migrationSpeed << endl; 
		return; 
	}

	/* --------------------------- PROLIFERATION RATE --------------------------- */
	/* NOTE: this function will not be used in IVDBM-ABM (stem cell version)      */
	void Agent::calculateProliferationRate(){
		/* Change in cell population (% of initial population) for gel with given %Alg(w/w) at time t (hours)
		 * 		   Cell population = (-a*Alg(w/w) + b)*ln(t_hours) + (c * Alg(w/w) + d)
		 *
		 *         Finite cell lines undergo linear then logarithmic cell growth. 
		 *         High Alg content and elastic moudlus hydrogels are favorable to cell attachment and proliferation  
		 */
					
		// Calculate population as percent of initial population at current (tick = t) and previous call (tick =t-1)
		float t = (WHWorld::reportHour());              // hours elapsed at tick t
		float tMinusOne = (WHWorld::clock-1)*0.5;       // hours elapsed at tick t-1
			
		#ifdef CALIBRATION
			float cellPopulation_t = 		 Agent::CaAlgProlif[0] - Agent::CaAlgProlif[1]*(Agent::agentWorldPtr->Alg_wv) - Agent::CaAlgProlif[2]*(Agent::agentWorldPtr->pXL) + Agent::CaAlgProlif[3]*(t)*(Agent::agentWorldPtr->Alg_wv) 		 + Agent::CaAlgProlif[4]*(Agent::agentWorldPtr->Alg_wv)*(Agent::agentWorldPtr->pXL); 
			float cellPopulation_tMinusOne = Agent::CaAlgProlif[0] - Agent::CaAlgProlif[1]*(Agent::agentWorldPtr->Alg_wv) - Agent::CaAlgProlif[2]*(Agent::agentWorldPtr->pXL) + Agent::CaAlgProlif[3]*(tMinusOne)*(Agent::agentWorldPtr->Alg_wv) + Agent::CaAlgProlif[4]*(Agent::agentWorldPtr->Alg_wv)*(Agent::agentWorldPtr->pXL); 

		#else 
			float cellPopulation_t =13.06 - 3.69*(Agent::agentWorldPtr->Alg_wv) - 101.18*(Agent::agentWorldPtr->pXL) + 0.056*(t)*(Agent::agentWorldPtr->Alg_wv) + 30.44*(Agent::agentWorldPtr->Alg_wv)*(Agent::agentWorldPtr->pXL); //14.635 + 6.95*(t) - 1.04*(t)*(Agent::agentWorldPtr->pXL) - 2.65*(Agent::agentWorldPtr->Alg_wv)*(Agent::agentWorldPtr->pXL); 
			float cellPopulation_tMinusOne = 13.06 - 3.69*(Agent::agentWorldPtr->Alg_wv) - 101.18*(Agent::agentWorldPtr->pXL) + 0.056*(tMinusOne)*(Agent::agentWorldPtr->Alg_wv) + 30.44*(Agent::agentWorldPtr->Alg_wv)*(Agent::agentWorldPtr->pXL); //14.635 + 6.95*(tMinusOne) - 1.04*(tMinusOne)*(Agent::agentWorldPtr->pXL) - 2.65*(Agent::agentWorldPtr->Alg_wv)*(Agent::agentWorldPtr->pXL);

		#endif

		// Change in cell population between current and previous time point, as percent of initial population:
		float deltaCellPopulation = (cellPopulation_t < 0 || cellPopulation_t - cellPopulation_tMinusOne < 0)? 0: cellPopulation_t - cellPopulation_tMinusOne; 

		// Calculate % of current chondrocytes that need to proliferate during current tick:
		Agent::proliferationRate = 100*(deltaCellPopulation*Agent::agentWorldPtr->initialCells[0])/(Agent::agentWorldPtr->chonds.size()); 
		cout << " Proliferation Rate (%) = " << Agent::proliferationRate << endl;         
		return; 
	}

	/* ----------------------------- VIABILITY RATE ----------------------------- */
	void Agent::calculateViabilityRate(){
		/* Calculate viability Rate (%) given current time (day)
		 * 		  Viability Rate (%) = a * ln(t_day) + b
		 *
		 *        Ca-Alg crosslinked cells are non-toxic and maintain cell viability and phenotype over time       
		 *        Ratio of Live to dead cells remains relatively constant over time       
		 */

		#ifdef CALIBRATION
			float viability = Agent::CaAlgViability[0]*log(WHWorld::reportDay()) + Agent::CaAlgViability[1];		
		#else
			float viability = 0.7321*log(WHWorld::reportDay()) + 97.452; 
		#endif

		if (viability < 0) viability = 0; 
		else if (viability > 100) viability = 100; 
		Agent::viabilityRate = viability; 

		cout << " Viability Rate (%)= " << Agent::viabilityRate << endl; 
		return; 
	}

	/* --------------------------- ECM SYNTHESIS RATE --------------------------- */
	void Agent::calculateECMSynthesisRate(){
		/* Cellular Activity and ECM synthesis dependent on Elastic Modulus (E) and pore size (poreWidth)
		 * 		  Collagen per cell = a*(b*E + c) + d*(-e*poreWidth + f)
		 * 		  Aggrecan per cell = g*(h*E + i) + j*(-k*poreWidth + l)
		 * 		  sGAG per cell = m*(n*E + o) + p*(-q*poreWidth + r)
		 *        
		 *        sGAG (HA) production increased with Elastic modulus 
		 *        Aggrecan production decreased with Elastic modulus 
		 *        Collagen production depends on pore size and not Elastic modulus 
		 */

		#ifdef CALIBRATION
			Agent::collagenSynthRate = Agent::CollagenSynth[0]*(Agent::CollagenSynth[1]*WHWorld::reportDay() + Agent::CollagenSynth[2]);
			Agent::aggrecanSynthRate = Agent::AggrecanSynth[0]*(Agent::AggrecanSynth[1]*WHWorld::reportDay() + Agent::AggrecanSynth[2]);
		#else		
			Agent::collagenSynthRate = 10*(6.45*WHWorld::reportDay() + 3.6);//10*(6.45*WHWorld::reportDay() + 3.6); // Agent::collagenSynthRate = 8922 - 255.4*mesh - 0.02429*Agent::agentWorldPtr->E; //Metabolism of the extracellular matrix formed by intervertebral disc cells cultured in alginate 
			Agent::aggrecanSynthRate = 20*(38*WHWorld::reportDay() + 16.6);//10*(38*WHWorld::reportDay() + 16.6);  // Agent::aggrecanSynthRate = 183.4 - 4.702*mesh - 0.0009759*Agent::agentWorldPtr->E;//metabolism of the extracellular matrix formed by intervertebral disc cells cultured in alginate		
		#endif
		cout << "        Collagen Synthesis Rate = " << Agent::collagenSynthRate << endl; 
		cout << "        Agggrecan Synthesis Rate = " << Agent::aggrecanSynthRate << endl; 
		return; 
	}

	void Agent::chondrocyteCaAlgBehavior(){
		if (WHWorld::clock == 0){
			Agent::calculateMigrationSpeed(); 
			Agent::calculateECMSynthesisRate();

		} else if (WHWorld::clock > 0) {
			Agent::calculateProliferationRate();
			Agent::calculateViabilityRate();
		}
	return; 
	}
#endif // MODEL_SCAFFOLD

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
	
	int tid = 0;
	#ifdef _OMP
		tid = omp_get_thread_num();	 // Get thread id in order to access the seed that belongs to this thread
	#endif

	int trial = 0;
	do {
		// Pick a neighbor to move to at random:
		#ifndef MODEL_3D
			int i = rand_r(&(agentWorldPtr->seeds[tid])) % 8;
		#else
			int i = rand_r(&(agentWorldPtr->seeds[tid])) % 27;
		#endif
			int dx = Agent::dX[neighbor[i]];
			int dy = Agent::dY[neighbor[i]];
			int dz = Agent::dZ[neighbor[i]];

		// Calculate target index:
		int newindex = (x + dx) + (y + dy)*nx + (z + dz)*nx*ny;

		// If the z-direction movement is invalid, pick a new neighbor
		if (z + dz < 0 || z + dz >= nz) continue;

		// If trying to move off the side boundaries, die:
		if (x + dx < 0 || x + dx >= nx || y + dy < 0 || y + dy >= ny) {
			this->life[write_t] = 0;
			this->die();
			return;
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
			int randInt = rand_r(&(agentWorldPtr->seeds[tid]))%(xtarget.size());
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

float Agent::meanNeighborChem(int chemIndex) {
	int totalchemical = 0, numberofpatches = 0;
  	
	// Location of agent in x,y,z dimensions of world.
	int x = this->ix[read_t];
	int y = this->iy[read_t];
	int z = this->iz[read_t];
  	
	// Number of patches in x,y,z dimensions of world
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	// Count number of chemicals of type chemIndex in all neighbors inside world dimensions:
	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
				if (x + dX < 0 || x + dX >= nx || y + dY < 0 || y + dY >= ny || z + dZ < 0 || z + dZ >= nz) continue;
				int in = (x + dX) + (y + dY)*nx + (z + dZ)*nx*ny;
				totalchemical += Agent::agentWorldPtr->chemAllocation[chemIndex][in];
				numberofpatches++;
			}
		}
	}
	return totalchemical/numberofpatches;
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

bool Agent::moveToHighestChem(int chemIndex) {
	int read_index;
	
	// Check if the location has been modified in this tick
	if (isModified(this->index)) read_index = write_t;	// If it has, work off of the intermediate value
	else read_index = read_t;							// If it has NOT, work off of the original value

    // Location of agent in x,y,z dimensions of world.
	int ix = this->ix[read_index];
	int iy = this->iy[read_index];
	int iz = this->iz[read_index];
	int index = this->index[read_index];
    
	// Number of patches in x,y,z dimensions of world.
	int nx = Agent::nx;
	int ny = Agent::ny;
	int nz = Agent::nz;

	double highestchem = Agent::agentWorldPtr->chemAllocation[chemIndex][index];
	int dx = 0, dy = 0, dz = 0;

	#ifdef MODEL_SCAFFOLD
        int radius = Agent::migrationSpeed; // Cells in Ca-Alg can move up to 'migrationSpeed' patches radius per tick 

        // Find the neighbor inside world dimensions with the highest concentration of chemical of type chemIndex:
		for (int deltaz = -radius; deltaz <= radius; deltaz++) {
			for (int deltay = -radius; deltay <= radius; deltay++) {
				for (int deltax = -radius; deltax <= radius; deltax++) {
					if (ix + deltax < 0 || ix + deltax >= nx || iy + deltay < 0 || iy + deltay >= ny || iz + deltaz < 0 || iz + deltaz >= nz) continue;
					
					int in = (ix + deltax) + (iy + deltay)*nx + (iz + deltaz)*nx*ny;
					if (Agent::agentWorldPtr->chemAllocation[chemIndex][in] > highestchem) {
						highestchem = Agent::agentWorldPtr->chemAllocation[chemIndex][in];
						dx = deltax;
						dy = deltay;
						dz = deltaz;
					}
				}
			}
		}

	#else
		// Find the neighbor inside world dimensions with the highest concentration of chemical of type chemIndex:
		for (int deltaz = -1; deltaz <= 1; deltaz++) {
			for (int deltay = -1; deltay <= 1; deltay++) {
				for (int deltax = -1; deltax <= 1; deltax++) {
					if (ix + deltax < 0 || ix + deltax >= nx || iy + deltay < 0 || iy + deltay >= ny || iz + deltaz < 0 || iz + deltaz >= nz) continue;
					
					int in = (ix + deltax) + (iy + deltay)*nx + (iz + deltaz)*nx*ny;
					if (Agent::agentWorldPtr->chemAllocation[chemIndex][in] > highestchem) {
						highestchem = Agent::agentWorldPtr->chemAllocation[chemIndex][in];
						dx = deltax;
						dy = deltay;
						dz = deltaz;
					}
				}
			}
		}
	#endif

	// Move to the neighbor with the highest concentration of chemical of type chemIndex if it is not the agent's current patch:
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
}