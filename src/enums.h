/* 
 * File: enums.h
 * 
 * File Contents: Contains enum declarations
 *
 * Author: Alireza Najafi-Yazdi
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#ifndef ENUMS_H
#define	ENUMS_H
#define REAL double

// Types of agents & ECM managers:
enum agent_t { 
  unoccupied = -1,      
  cell = 1,
  stem = 2,
  progen = 3,
  np = 4,
  // achondrocyte = 6,     
  orig_coll = 9,    // Original collagen
  new_coll = 10,    // New collagen
  frag_coll = 11,   // Fragmented collagen
  orig_agg = 12,    // Original aggrecan
  new_agg = 13,     // New aggrecan
  frag_agg = 14,    // Fragmented aggrecan
  oha = 15,  
  nha = 16,  
  fha = 17  
};

/* Colors of patches and agents for visualization with Paraview 3.0 */
enum color_t {
  ccell = 105,
  cstem = 100,
  cprogen = 90,
  cnp = 80,

  ccollagen = 139,
  caggrecan = 14,  
  
  cfcollagen = 0,
  cfaggrecan = 0,

  cHA = 0,
  cfHA = 0, 
  cnothing = 0,          
  cdamage = 0,
  cunidentifiable = 0, 

  cCaAlg = 119 
};

/** Species ids and matching grid channel indices (see ChemicalEnvironment). */
enum chemical_t {
  /** Species identifiers for agent/world API. */
  TNF = 0,
  TGF = 1,
  IL1beta = 2,

  pTNF = 0,
  pTGF = 1,
  pIL1beta = 2,

  dTNF = 3,
  dTGF = 4,
  dIL1beta = 5,
  pcellgrad = 6,
};

// Types of patches:
enum patches_t {
  nothing = 0,          
  damage = 4,
  unidentifiable = 5,
  CaAlg = 10
};

// Time points within a tick
enum readwrite_t {
  read_t = 0,  // Start of a tick
  write_t = 1  // End of a tick
}; 

//enum peptides {
//	MAL,
//	CHAD,
//	hA5G26,
//	IKVAV
//};

#endif	/* ENUMS_H */