/* 
 * File: ECM.cpp
 *
 * File Contents: Contains ECM class
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#include "ECM.h"
#include "../World/Usr_World/biomaterialWorld.h"
#include "../enums.h"
#include <iostream>
#include <vector>
#include <cstdlib>                                      
#include <stdio.h>                                     
#include <string.h>                                     
#include <algorithm>

using namespace std;

//FIXME: Update max num of ECM on each patch 
Patch* ECM::ECMPatchPtr = NULL; 
BMWorld* ECM::ECMWorldPtr = NULL;
int ECM::maxcollagen = 620*10^9;  
int ECM::maxaggrecan = 500*10^9;  
int ECM::maxHA = 0;
int ECM::dx[27] = {-1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1};
int ECM::dy[27] = {-1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1};
int ECM::dz[27] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1};
int ECM::d[27] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};    

ECM::ECM() {}

ECM::ECM(int x, int y, int z, int index) {
	this->dirty = true;
	this->request_dirty = true;
	this->dirty_from_neighbors = true;
	this->indice[0] = x;
	this->indice[1] = y;
	this->indice[2] = z;
	this->index = index;
	this->rng = abm::rng::Stream(static_cast<uint32_t>(index), 0);
	this->empty[read_t] = true;
	this->ocollagen[read_t] = 0;
	this->ncollagen[read_t] = 0;
	this->fcollagen[read_t] = 0;
	this->oaggrecan[read_t] = 0;
	this->naggrecan[read_t] = 0;
	this->faggrecan[read_t] = 0;
	this->HA[read_t] = 0;
	this->fHA[read_t] = 0;
	this->fcollDangerSignal[read_t] = false;
	this->faggDangerSignal[read_t] = false;
	this->fHADangerSignal[read_t] = false;
	this->scarIndex[read_t] = false;
	this->empty[write_t] = true;
	this->ocollagen[write_t] = 0;
	this->ncollagen[write_t] = 0;
	this->fcollagen[write_t] = 0;
	this->oaggrecan[write_t] = 0;
	this->naggrecan[write_t] = 0;
	this->faggrecan[write_t] = 0;
	this->HA[write_t] = 0;
	this->fHA[write_t] = 0;
	this->fcollDangerSignal[write_t] = false;
	this->faggDangerSignal[write_t] = false;
	this->fHADangerSignal[write_t] = false;
	this->scarIndex[write_t] = false;

	memset(requestfcollagen, 0, 27*sizeof(int));
	memset(requestfaggrecan, 0, 27*sizeof(int));
	memset(requestfHA, 0, 27*sizeof(int));
}

ECM::~ECM() {}

void ECM::set_dirty() {
	#ifdef OPT_ECM
		this->dirty = true;
	#endif
}

void ECM::set_request_dirty() {
	#ifdef OPT_ECM
		this->request_dirty = true;
	#endif
}

void ECM::set_dirty_from_neighbors() {
	#ifdef OPT_ECM
		this->dirty_from_neighbors = true;
	#endif
}

void ECM::reset_dirty() {
	#ifdef OPT_ECM
		this->dirty = false;
	#endif
}

void ECM::reset_request_dirty() {
	#ifdef OPT_ECM
		this->request_dirty = false;
	#endif
}

void ECM::reset_dirty_from_neighbors() {
	#ifdef OPT_ECM
		this->dirty_from_neighbors = false;
	#endif
}

void ECM::decrement(int n) {
	--n;
}

void ECM::ECMFunction() {

	/* -------------------------------------------------------------------------- */
	/*                                DAMAGE REPAIR                               */
	/* -------------------------------------------------------------------------- */
	// New ECM proteins can repair damage on their patch:
	if (this->ncollagen[read_t] > 0 || this->naggrecan[read_t] > 0) this->repairDamage();

	/* -------------------------------------------------------------------------- */
	/*                                    DEATH                                   */
	/* -------------------------------------------------------------------------- */
   
	/* Each hyaluronan life decreases at each time step and they die naturally
	 * Collagen & aggrecan have 'infinite life' because their half lives are on the order of years. */
	#ifdef OPT_ECM
		for_each(HAlife.begin(), HAlife.end(), decrement);
	#else
		for (int i = 0; i < HAlife[read_t].size(); i++) {
		HAlife[write_t][i] = HAlife[read_t][i] - 1;
		}
	#endif

	/* -------------------------------------------------------------------------- */
	/*                              DANGER SIGNALLING                             */
	/* -------------------------------------------------------------------------- */
	// Fragmented ECM proteins can signal danger one time:
	if (fcollDangerSignal[read_t] == true || faggDangerSignal[read_t] == true) {
		this->set_dirty();
		this->dangerSignal();
		fcollDangerSignal[write_t] = false;
		faggDangerSignal[write_t] = false;
	}

	/* -------------------------------------------------------------------------- */
	/*                               SCAR FORMATION                               */
	/* -------------------------------------------------------------------------- */
	// Original collagen can create a scar if above threshold of 100:
	if (ncollagen[read_t] >= 10) {  // TODO after sensitivity
		this->set_dirty();
		scarIndex[write_t] = true;                                 // FIXME
	}
}

void ECM::repairDamage() {
  	// Location of neighbor in x,y,z dimensions of the world:
	int tempX, tempY, tempZ, tempIndex;
	int nx, ny, nz;
	Patch* tempPatchPtr;

  	// Number of patches in x,y,z dimensions of the world:
	nx = ECM::ECMWorldPtr->nx;
	ny = ECM::ECMWorldPtr->ny;
	nz = ECM::ECMWorldPtr->nz;

	//Repair damage on current and Neighbor Patches
	for (int dZ = -1; dZ <= 1; dZ++) {
		for (int dY = -1; dY <= 1; dY++) {
			for (int dX = -1; dX <= 1; dX++) {
        		// Location of patch to be repaired in x,y,z dimensions of world.
				tempX = this->indice[0] + dX;
				tempY = this->indice[1] + dY;
				tempZ = this->indice[2] + dZ;

        		// Try a new patch if this one is outside the world dimensions.
				if (tempX < 0 || tempX >= nx || tempY < 0 || tempY >= ny || tempZ < 0 || tempZ >= nz) continue;

        		// Get access to the patch on which this ECM manager resides
				tempIndex = tempX + tempY*nx + tempZ*(nx)*(ny);
				tempPatchPtr = &(ECM::ECMPatchPtr[tempIndex]);

				#ifdef MODEL_SCAFFOLD
					// Repair Damage, but do not replace CaAlg scaffold 
					if (tempPatchPtr->damage[write_t] > 0) {

						// Data race allowed, since we're just overwriting values
						tempPatchPtr->dirty = true;
						tempPatchPtr->damage[write_t] = 0;
						tempPatchPtr->health[write_t] = 100;                                

						if (tempPatchPtr->type[read_t] != CaAlg){
							tempPatchPtr->type[write_t] = nothing; 
							tempPatchPtr->color[write_t] = cnothing; 
						}
					}
				#else
					// Repair Damage
					if (tempPatchPtr->damage[write_t] > 0) {
						// Data race allowed, since we're just overwriting values:
						tempPatchPtr->dirty = true;
						tempPatchPtr->damage[write_t] = 0;
						tempPatchPtr->type[write_t] = tissue;
						tempPatchPtr->color[write_t] = ctissue;
						tempPatchPtr->health[write_t] = 100;
						// cout << " repair damage at " <<tempIndex << endl; 
					}
				#endif
			}
		}
	}
}

void ECM::dangerSignal() {
	this->ECMPatchPtr[this->index].dirty = true;
	this->ECMPatchPtr[this->index].damage[write_t]++;
	this->ECMPatchPtr[this->index].health[write_t] = 0;
	#ifndef MODEL_SCAFFOLD
		this->ECMPatchPtr[this->index].color[write_t] = cdamage;
	#endif

	// Activated cells can remove newly created damage if they are present.
	if (ECMPatchPtr[this->index].isOccupied() == false) return;
}

void ECM::fragmentNCollagen() {
  	// Distance to neighbor in x,y,z dimensions of the world:
	int dX, dY, dZ;
	int newfragments = 0;

  	// Location of ECM manager in x,y,z dimensions of the world:
	int ix = indice[0];
	int iy = indice[1];
	int iz = indice[2];
  
  	// Number of patches in x,y,z dimensions of the world:
	int nx = ECM::ECMWorldPtr->nx;
	int ny = ECM::ECMWorldPtr->ny;
	int nz = ECM::ECMWorldPtr->nz;

  	// Alert change in status of original collagen
	if (this->ncollagen[write_t] > 0) {
		this->set_dirty();
		this->set_request_dirty();
	}

	// Request hatching two fragmented collagens on neighboring patches for each original collagen:
	while (this->ncollagen[write_t] > 0) {
		this->ncollagen[write_t]--;
		newfragments = 0;

		// Request one fragmented collagen at a random inbounds neighboring patch. //TODO(Caroline) Might want to make the radius = 2
		int d[27];
		std::copy(&ECM::d[0], &ECM::d[27], &d[0]);
		abm::rng::shuffle(this->rng, &d[0], &d[27]);
		for (int i = 0; i < 27; i++) {
			dX = dx[d[i]];
			dY = dy[d[i]];
			dZ = dz[d[i]];
			if (newfragments >= 2) break;
			if (ix + dX < 0 || ix + dX >= nx || iy + dY < 0 || iy + dY >= ny || iz + dZ < 0 || iz + dZ >= nz) continue;  //'dX + dY*3 + dZ*3*3 + 13' determines which neighbor
			this->requestfcollagen[dX + dY*3 + dZ*3*3 + 13]++;
			int in = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
     		
			// Alert change in status of collagen on this patch
			this->ECMWorldPtr->worldECM[in].set_dirty_from_neighbors();
			newfragments++;
		}
	}
	this->isEmpty();
}

void ECM::fragmentNAggrecan() {
	// Distance to neighbor in x,y,z dimensions of the world:
	int dX, dY, dZ;
	int newfragments = 0;
	
	// Location of ECM manager in x,y,z dimensions of the world:
	int ix = indice[0];
	int iy = indice[1];
	int iz = indice[2];
	
	// Number of patches in x,y,z dimensions of the world:
	int nx = ECM::ECMWorldPtr->nx;
	int ny = ECM::ECMWorldPtr->ny;
	int nz = ECM::ECMWorldPtr->nz;

  	// Alert change in status of original aggrecan
	if (this->naggrecan[write_t] > 0) {
		this->set_dirty();
		this->set_request_dirty();
	}

	// Request hatching two fragmented aggrecans on neighboring patches for each original aggrecan:
	while (this->naggrecan[write_t] > 0) {
		this->naggrecan[write_t]--;
		newfragments = 0;

		// Request one fragmented aggrecan at a random inbounds neighboring patch. // TODO(Caroline) Might want to make the radius = 2
		int d[27];
		std::copy(&ECM::d[0], &ECM::d[27], &d[0]);
		abm::rng::shuffle(this->rng, &d[0], &d[27]);
		for (int i = 0; i < 27; i++) {
			dX = dx[d[i]];
			dY = dy[d[i]];
			dZ = dz[d[i]];
			if (newfragments >= 2) break;
			if (ix + dX < 0 || ix + dX >= nx || iy + dY < 0 || iy + dY >= ny || iz + dZ < 0 || iz + dZ >= nz) continue; //'dX + dY*3 + dZ*3*3 + 13' determines which neighbor
			this->requestfaggrecan[dX + dY*3 + dZ*3*3 + 13]++;
			int in = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
      	
			// Alert change in status of aggrecan on this patch
			this->ECMWorldPtr->worldECM[in].set_dirty_from_neighbors();
			newfragments++;
		}
	}
	this->isEmpty();
}

void ECM::fragmentHA() {
  	// Distance to neighbor in x,y,z dimensions of the world:
	int dX, dY, dZ;
	int newfragments = 0;
  	
	// Location of ECM manager in x,y,z dimensions of the world:
	int ix = indice[0];
	int iy = indice[1];
	int iz = indice[2];
  	
	// Number of patches in x,y,z dimensions of the world:
	int nx = ECM::ECMWorldPtr->nx;
	int ny = ECM::ECMWorldPtr->ny;
	int nz = ECM::ECMWorldPtr->nz;

  	// Alert change in status of hyaluronan
	if (HA[write_t] > 0) {
		this->set_dirty();
		this->set_request_dirty();
	}

	// Request hatching two fragmented hyaluronans on neighboring patches for each hyaluronan :
	while (HA[write_t] > 0) {
		this->HA[write_t]--;
		newfragments = 0;

   // Kills hyaluronan that is to be fragmented:
	#ifdef OPT_ECM
		int randomHA = abm::rng::uniform_int(this->rng, HAlife.size());
		HAlife[randomHA] = 0;
	#else
		int randomHA = abm::rng::uniform_int(this->rng, HAlife[write_t].size());
		HAlife[write_t][randomHA] = 0;
	#endif

	// Request one fragmented hyaluronan at a random inbounds neighboring patch. TODO(Caroline) Might want to make the radius = 2
		int d[27];
		std::copy(&ECM::d[0], &ECM::d[27], &d[0]);
		abm::rng::shuffle(this->rng, &d[0], &d[27]);
		for (int i = 0; i < 27; i++) {
			dX = dx[d[i]];
			dY = dy[d[i]];
			dZ = dz[d[i]];
			if (newfragments >= 2) break;
			if (ix + dX < 0 || ix + dX >= nx || iy + dY < 0 || iy + dY >= ny || iz + dZ < 0 || iz + dZ >= nz) continue; //'dX + dY*3 + dZ*3*3 + 13' determines which neighbor in radius 1
			this->requestfHA[dX + dY*3 + dZ*3*3 + 13]++;
			int in = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
      	
			// Alert change in status of hyaluronan on this patch
			this->ECMWorldPtr->worldECM[in].set_dirty_from_neighbors();
			newfragments++;
		}
	}
	this->isEmpty();
}

void ECM::updateECM() {
	int in;   // Patch row major index of neighbor:
	int fcollagenrequest = 0, faggrecanrequest = 0, fHArequest = 0;   // Amount of requested fragmented ECM proteins:

  	// Location of ECM manager in x,y,z dimensions of the world:
	int ix = this->indice[0];
	int iy = this->indice[1];
	int iz = this->indice[2];
  
  	// Number of patches in x,y,z dimensions of the world:
	int nx = ECM::ECMWorldPtr->nx;
	int ny = ECM::ECMWorldPtr->ny;
	int nz = ECM::ECMWorldPtr->nz;

  /*************************************************************************
   * FRAGMENTED ECM REQUESTS                                               *
   *************************************************************************/
	// Iterate through neighboring patches, count any fcollage/aggrecan/HA requests for ECM manager
	#ifdef OPT_ECM
		if (this->dirty_from_neighbors) { 	// Only process requests if a neighbor has indicated that it's made a fragment request to this ECM manager
	#endif

#ifdef ECM_UNROLL_LOOP
    // Distance to neighbors in x,y,z dimensions of the world:
  	int dX[8] = {-1, 0, 1, -1, 1, -1, 0, 1}; 
    int dY[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    int dZ = 0;  // 2D
	int targetZ = iz + dZ;
    for (int i = 0; i < 8; i++) {
        int targetX = ix + dX[i];
		int targetY = iy + dY[i];
        // Only consider requests from neighbors that are inside the world.
		    if (!(targetX < 0 || targetX >= nx || targetY < 0 || targetY >= ny || targetZ < 0 || targetZ >= nz)) {
	    		in = targetX + targetY*nx + targetZ*nx*ny;
	    		ECM* neighborECMPtr = &(ECMWorldPtr->worldECM[in]);
          		
				/* self_neighbor_in is the index of this ECM manager in the list of its neighbor's list of neighbors. */
	    		int self_neighbor_in = (-dX[i]) + (-dY[i])*3 + (-dZ)*3*3 + 13;
    			fcollagenrequest += neighborECMPtr->requestfcollagen[self_neighbor_in];
    			faggrecanrequest += neighborECMPtr->requestfaggrecan[self_neighbor_in];
	    	}
    }
#else
    // dX, dY, dZ are distance to neighbors in x,y,z dimensions of the world
		for (int dX = -1; dX < 2; dX++) {
			for (int dY = -1; dY < 2; dY++) {
				for (int dZ = -1; dZ < 2; dZ++) {
          // Try another neighbor if this one is out of bounds
					if (ix + dX < 0 || ix + dX >= nx || iy + dY < 0 || iy + dY >= ny || iz + dZ < 0 || iz + dZ >= nz) continue;
					in = (ix + dX) + (iy + dY)*nx + (iz + dZ)*nx*ny;
					int a = ECMWorldPtr->worldECM[in].requestfcollagen[(-dX) + (-dY)*3 + (-dZ)*3*3 + 13];
					int b = ECMWorldPtr->worldECM[in].requestfaggrecan[(-dX) + (-dY)*3 + (-dZ)*3*3 + 13];
					int c = ECMWorldPtr->worldECM[in].requestfHA[(-dX) + (-dY)*3 + (-dZ)*3*3 + 13];
					if (a != 0) {
						//cout << "  fcollagen requested at " << ix + iy*nx + iz*nx*ny << " by " << in << endl;
					}
					if (b != 0) {
						//cout << "  faggrecan requested at " << ix + iy*nx + iz*nx*ny << " by " << in << endl;
					}
					fcollagenrequest += a;
					faggrecanrequest += b;
				}
			}
		}
#endif

		// Fragmented ECM requests can only be accepted if there is enough space. // TODO(Kim): INSERT REF?
		if (ocollagen[read_t] + ncollagen[read_t] + fcollagen[read_t] + fcollagenrequest > maxcollagen) {
			cout << "  Error fcollagen request" << endl;
		} else if (oaggrecan[read_t] + naggrecan[read_t] + faggrecan[read_t] + faggrecanrequest > maxaggrecan) {
			cout << "  Error faggrecan request" << endl;
		} else {
      		// Fragmented ECM proteins serve as danger signals once
			this->fcollagen[write_t] += fcollagenrequest;
			this->fcollDangerSignal[write_t] += fcollagenrequest;
			this->faggrecan[write_t] += faggrecanrequest;
			this->faggDangerSignal[write_t] += faggrecanrequest;
			this->isEmpty();
		}
		this->fcollagen[read_t] = this->fcollagen[write_t];
		this->faggrecan[read_t] = this->faggrecan[write_t];

#ifdef OPT_ECM
	}	// if (this->dirty_from_neighbors)

  /*************************************************************************
   * READ/WRITE SYNCHRONIZATION                                            *
   *************************************************************************/
	// Only synchronize the read and write entries if there's a change in value
	if (this->dirty) {
#endif
		this->empty[read_t] = this->empty[write_t];
		this->ocollagen[read_t] = this->ocollagen[write_t];
		this->ncollagen[read_t] = this->ncollagen[write_t];
		this->oaggrecan[read_t] = this->oaggrecan[write_t];
		this->naggrecan[read_t] = this->naggrecan[write_t];
		this->HA[read_t] = this->HA[write_t];

		// Convert ocollagen (tropocollagen monomer) to ncollagen (polymer)
		while (this->ocollagen[read_t] > 1) {
			this->ocollagen[read_t] -= 2;
			this->ocollagen[write_t] -= 2;
			this->ncollagen[read_t] += 1;
			this->ncollagen[write_t] += 1;
		}

		while (this->oaggrecan[read_t] > 1) {
			this->oaggrecan[read_t] -= 2;
			this->oaggrecan[write_t] -= 2;
			this->naggrecan[read_t] += 1;
			this->naggrecan[write_t] += 1;
		}

    // Get rid of all dead hyaluronans
	#ifdef OPT_ECM
		vector<int> *vec = &(this->HAlife);
	#else
		vector<int> *vec = &(this->HAlife[write_t]);
	#endif
		vector<int>::iterator first	= vec->begin();
		for (vector<int>::iterator it = first; it != vec->end();) {
			int life = *it;
			if (life <= 0) it = vec->erase(it);
			else ++it;	
		}

		this->fcollDangerSignal[read_t] = this->fcollDangerSignal[write_t];
		this->faggDangerSignal[read_t] = this->faggDangerSignal[write_t];
		this->scarIndex[read_t] = this->scarIndex[write_t];
	#ifdef OPT_ECM
		}	// if (this->dirty)
	#endif

	// If number of oECM exceed threshold, replace with nECM define ocollagen to ncollagen (collagen molecule) conversion threshold, after sensitivity analysis 

	// Remove all dirty flags
	this->reset_dirty();
	this->reset_dirty_from_neighbors();
}

void ECM::isEmpty() {
	int totalcollagen = ocollagen[write_t] + ncollagen[write_t] + fcollagen[write_t] ;
	int totalaggrecan = oaggrecan[write_t] + naggrecan[write_t] + faggrecan[write_t];
		if (totalcollagen + totalaggrecan == 0) { 
		this->empty[write_t] = true;
	} else {
		this->empty[write_t] = false;
	}
}

void ECM::resetrequests() {
#ifdef OPT_ECM
	// Only clear the requests if there are any in this tick
	if (this->request_dirty) {
#endif
		memset(this->requestfcollagen, 0, 27*sizeof(int));
		memset(this->requestfaggrecan, 0, 27*sizeof(int));
#ifdef OPT_ECM
	}
#endif

	this->reset_request_dirty();
}