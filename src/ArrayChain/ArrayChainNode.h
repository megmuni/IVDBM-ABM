/*
 * ArrayChainNode.h
 *
 * File Contents: Contains the ArrayChainNode class
 *
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#ifndef ARRAYCHAINNODE_H_
#define ARRAYCHAINNODE_H_

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <vector>
#include "../common.h"

class Chondrocyte;

/*
 * ARRAYCHAINNODE CLASS DESCRIPTION:     The ArrayChainNode class manages agent/cell location information.
 *                                       It is used to keep track of free elements that agents could occupy, elements that
 *                                       agents have already occupied, and to update this information efficiently
 *                                       with multiple threads (including allocating more space for new data). 
 */
template <typename T>
class ArrayChainNode {
 public:
	typedef void (*DestroyElementFunc) (T& val);
  /*
   * Description: Default ArrayChainNode constructor. 
   *
   * Return: void
   *
   * Parameters: void
   */
	ArrayChainNode() {
		data = NULL;
		numFreeApprox = 0;
		capacity = 0;
		next = NULL;
		destroyFunc = NULL;
	}

  /*
   * Description: ArrayChainNode constructor. 
   *
   * Return: void
   *
   * Parameters: capacity  -- Maximum number of elements allowed in this node
   *             desFunc   -- Function pointer for destroying an entry in the array
   */
	ArrayChainNode(sizeType capacity, DestroyElementFunc desFunc) {
		this->numFreeApprox = capacity;
		this->capacity = capacity;
		allocateData(capacity);		// Allocate memory for data
		this->next = NULL; 		    // Set pointer to next node to NULL
		destroyFunc = desFunc;		// Set function pointer to destroy element
	}

  /*
   * Description: Copy constructor for ArrayChainNode. 
   *
   * Return: void
   *
   * Parameters: source  -- Reference to ArrayChainNode to be copied
   */
	ArrayChainNode(const ArrayChainNode& source) {
		data = source.data;
		numFreeApprox = source.numFreeApprox;
		capacity = source.capacity;
		next = source.next;
		destroyFunc = source.destroyFunc;
	}

  /*
   * Description: Overloaded assigment operator for ArrayChainNode.
   *
   * Return: ArrayChainNode that has been copied over
   *
   * Parameters: source  -- Reference to ArrayChainNode to be copied
   */
	ArrayChainNode& operator=(const ArrayChainNode& source) {
		if (this == &source) return *this;
		
		data = source.data;
		numFreeApprox = source.numFreeApprox;
		capacity = source.capacity;
		next = source.next;
		destroyFunc = source.destroyFunc;
		return *this;
	}

  /*
   * Description: Empties an ArrayChainNode
   *
   * Return: void
   *
   * Parameters: void
   */
	void emptySelf() {
		if (destroyFunc) {
			for (int i = 0; i < capacity; i++) destroyFunc(data[i]);	
		} else {
			std::cerr << "		ArrayChainNode::emptySelf(): NULL destroyFunc" << std::endl;
		}
		free(data);
	}

  /*
   * Description: ArrayChainNode destructor. 
   *
   * Return: void
   *
   * Parameters: void
   */
	~ArrayChainNode() {
		emptySelf();
	}

  // Pointer to the next node in the chain
	ArrayChainNode<T>* next;

	/* -------------------------------------------------------------------------- */
	/*                                   Setters                                  */
	/* -------------------------------------------------------------------------- */
  /*
   * Description: Adds data to an existing ArrayChainNode at a given thread.
   *              NOTE: Not thread safe.
   *
   * Return: void
   *
   * Parameters: val  -- Data to be added (cell pointer)
   *             tid  -- Thread identification number
   */
	bool addData(T val, int tid) {
		// Check for locally known available element and add data there:
		if (freeLists[tid].size()) {
			sizeType index = freeLists[tid].back();  // read index of empty element
			freeLists[tid].pop_back();  // Remove element from the list of free ones
			data[index] = val;  // Add the data to the empty element
			return true;
		} else if (numFreeApprox) {  // No known spots, but could be some free
			data[capacity - numFreeApprox] = val;  // Add data to space at the end
			numFreeApprox--;  // Update number of free elements at the end
			return true;
		}
		// Return false for the chain to go to next node or allocate a new node
		return false;
	}

  /*
   * Description: Removes data from an existing ArrayChainNode at a given thread
   *
   * Return: void
   *
   * Parameters: nodeIndex  -- Index of the element where deletion will occur
   *             tid        -- Thread identification number
   */
	bool deleteData(sizeType nodeIndex, int tid) {
		if (nodeIndex > (capacity-1)) {
			return false;
		}
		// Delete element if user specified a destroy function
		if (destroyFunc) {
			destroyFunc(data[nodeIndex]);
		}
		data[nodeIndex] = NULL;
		addToFreeList(tid, nodeIndex);  // Add node to list of available elements
		return true;
	}

  /*
   * Description: Adds data to an existing ArrayChainNode at a given thread.
   *              NOTE: Thread safe.
   *
   * Return: void
   *
   * Parameters: val  -- Data to be added (cell pointer)
   *             tid  -- Thread identification number
   */
	bool concAddData(T val, int tid) {
		std::cout << "concAddData() not implemented yet" << std::endl;
		return false;
	}

	/* -------------------------------------------------------------------------- */
	/*                                   Getters                                  */
	/* -------------------------------------------------------------------------- */
  /*
   * Description: Gets approximate number of free elements in this ArrayChainNode
   *
   * Return: approximate number of free elements in this ArrayChainNode
   *
   * Parameters: void
   */
	sizeType getNumFreeApprox() {
		return numFreeApprox;
	}

  /*
   * Description: Gets the capacity of this ArrayChainNode
   *
   * Return: capacity of this ArrayChainNode
   *
   * Parameters: void
   */
	sizeType getCapacity() {
		return capacity;
	}

  /*
   * Description: Gets the data (cell pointer) at index
   *
   * Return: Data (cell pointer) at index
   *
   * Parameters: index  -- Index of the element to get data from
   */
	T& getDataAt(sizeType index) {
		return data[index];
	}

  /*
   * Description: Gets the next node in the ArrayChain
   *
   * Return: next node in ArrayChain
   *
   * Parameters: void
   */
	ArrayChainNode<T>* getNextNodePtr() {
		return next;
	}

  /*
   * Description: Gets actual size of the ArrayChainNode (occupied elements)
   *
   * Return: actual size of ArrayChainNode
   *
   * Parameters: void
   */
	sizeType actualSize() {
		sizeType sz = capacity - numFreeApprox;
		for (int i = 0; i < NUM_THREAD; i++) {
			sz -= freeLists[i].size();
		}
		return sz;
	}

 private:
  // Stores the data (cell pointers)
	T* data;
  
  /* Stores lists of indices of free elements in data.
   * Each vector contains best effort information local to each thread.
   * freeLists[tid] only contains indices of the elements in data that have
   * been emptied out by thread tid. Thus, each list or even all of them
   * combined may not be exhaustive. */
	std::vector<sizeType> freeLists[NUM_THREAD];

	sizeType numFreeApprox;   // Stores a number of contiguous elements at the end of the list. This is used when threads request for parallel batch insertion to data.
	sizeType capacity;        // Maximum number of elements allowed in this node
	DestroyElementFunc destroyFunc;   // Function pointer for destroying an entry in the array

	//  ArrayChainNode<T>* next; // Pointer to the next node in the chain
	//  void addToFreeList(int tid, sizeType index);  // index: what to add to list
 
	/* -------------------------------------------------------------------------- */
	/*                        Private Auxilliary Functions                        */
	/* -------------------------------------------------------------------------- */
  /*
   * Description: Allocates space for n T's (cell pointers)
   *
   * Return: void
   *
   * Parameters: n  -- Number of elements to make available for data
   */
	void allocateData(sizeType n) {
		data = (T*) malloc (n*sizeof(T));
		if (! data) {
			std::cout << "Error: Failed to allocate memory for " << n << "elements in allocateData()" << std::endl;
			exit(-1);
		}
		capacity = n;
	}

  /*
   * Description: Adds an element to the list of free elements for a given thread
   *
   * Return: void
   *
   * Parameters: tid    -- Thread identification number
   *             index  -- Index of the element to add to the free list
   */
	void addToFreeList(int tid, sizeType index) {
		freeLists[tid].push_back(index);
	}
};

#endif 