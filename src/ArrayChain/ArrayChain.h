/*
 * ArrayChain.h
 *
 * File Contents: Contains declarations for the ArrayChain class.
 *
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#ifndef ARRAYCHAIN_H_
#define ARRAYCHAIN_H_

#include "ArrayChainNode.h"

/*
 * ARRAYCHAIN CLASS DESCRIPTION:       The ArrayChain class manages agent/cell location information.
 *                                     An ArrayChain is composed of a chain of many objects of type ArrayChainNode.
 *                                     The class is used to add and remove data, create new ArrayChainNodes and print data.
 */
template<typename T>
class ArrayChain {
 public:
  typedef void (*PrintElementFunc) (T& val);
  typedef void (*DestroyElementFunc) (T& val);

  /*
   * Description: Default ArrayChain constructor. 
   *
   * Return: void
   *
   * Parameters: void
   */
  ArrayChain();

  /*
   * Description: ArrayChain constructor. 
   *
   * Return: void
   *
   * Parameters: nodeCap   -- Maximum number of elements allowed in an ArrayChainNode in this ArrayChain
   *             chainCap  -- Maximum number of elements in this chain
   *             pFunc     -- Function pointer for printing data in an element
   *             dFunc     -- Function pointer for destroying an element
   */
  ArrayChain(sizeType nodeCap, sizeType chainCap, PrintElementFunc pFunc, DestroyElementFunc dFunc);

  /*
   * Description: Copy constructor for ArrayChain. 
   *
   * Return: void
   *
   * Parameters: source  -- Reference to ArrayChain to be copied
   */
  ArrayChain(const ArrayChain& source);

  /*
   * Description: Overloaded assigment operator for ArrayChain.
   *
   * Return: ArrayChain that has been copied over
   *
   * Parameters: source  -- Reference to ArrayChain to be copied
   */
  ArrayChain& operator=(const ArrayChain& source);

  /*
   * Description: Virtual ArrayChain destructor. 
   *
   * Return: void
   *
   * Parameters: void
   */
  virtual ~ArrayChain();

	/* -------------------------------------------------------------------------- */
	/*                                   Setters                                  */
	/* -------------------------------------------------------------------------- */
  /*
   * Description: Adds data to the ArrayChain at a node for a given thread
   *              NOTE: Not thread safe.
   *
   * Return: void
   *
   * Parameters: val  -- Data to be added (cell pointer)
   *             tid  -- Thread identification number
   */
  bool addData(T val, int tid);

  /*
   * Description: Removes data from the ArrayChainNode element at globalIndex at a given thread.
   *
   * Return: void
   *
   * Parameters: globalIndex  -- Index of the element where deletion will occur
   *                             (Example: globalIndex=50, nodecapacity=20: tenth element of third node)
   *             tid          -- Thread identification number
   */
  bool deleteData(sizeType globalIndex, int tid);

  /*
   * Description: Adds data to the ArrayChain at a node for a given thread
   *              NOTE: Thread safe
   *
   * Return: void
   *
   * Parameters: val  -- Data to be added (cell pointer)
   *             tid  -- Thread identification number
   */
  bool concAddData(T val, int tid);

	/* -------------------------------------------------------------------------- */
	/*                                   Getters                                  */
	/* -------------------------------------------------------------------------- */
  /*
   * Description: Gets the data (cell pointer) from the element with globalIndex
   *
   * Return: Data (cell pointer) at globalIndex
   *
   * Parameters: globalIndex  -- Index of element to get data from
   *                             (Example: globalIndex=50, nodecapacity=20: tenth element of third node)
   */
  T& getDataAt(sizeType globalIndex);

  /*
   * Description: Gets the approximate size of the ArrayChain (number of elements - approximate free elements).
   *              Used to iterate through the ArrayChain.
   *
   * Return: Approximate size of ArrayChain (number of elements)
   *
   * Parameters: void
   */
  sizeType size();

  /*
   * Description: Gets actual size of the ArrayChainNode (occupied elements)
   *
   * Return: actual size of ArrayChainNode
   *
   * Parameters: void
   */
  sizeType actualSize();

	/* -------------------------------------------------------------------------- */
	/*                                   Printer                                  */
	/* -------------------------------------------------------------------------- */

  /*
   * Description: Prints the data at all elements in the ArrayChain.
   *              NOTE: To use this you must register a print function when the ArrayChain is initialized.
   *
   * Return: void
   *
   * Parameters: void
   */
  void print();

 private:
  ArrayChainNode<T>* head;   // Pointer to the first ArrayChainNode in the ArrayChain
  ArrayChainNode<T>* tail;   // Pointer to the last ArrayChainNode in the ArrayChain
  int numNodes;              // Number of nodes in this ArrayChain
  sizeType chainCapacity;    // Maximum number of elements in this ArrayChain
  sizeType nodeCapacity;     // Maximum number of elements allowed in an ArrayChainNode in this ArrayChain
  PrintElementFunc printFunc;       // Function pointer for printing data in an element of the ArrayChain
  DestroyElementFunc destroyFunc;   // Function pointer for destroying an element of the ArrayChain


	/* -------------------------------------------------------------------------- */
	/*                        Private Auxilliary Functions                        */
	/* -------------------------------------------------------------------------- */
  /*
   * Description: Converts a globalIndex (index of an element in the whole ArrayChain) to a localIndex (index of an element in its ArrayChainNode) 
   *              and updates the node number and the local index. 
   *              EXAMPLE: globalIndex=50, nodecapacity=20: tenth element of third node
   *
   * Return: void
   *
   * Parameters: globalIndex  -- Index of element in the ArrayChain
   *             nodeNum      -- Pointer to the node number that the element at globalIndex is in
   *             localIndex   -- Index of element in its ArrayChainNode
   */
  void globalToLocalIndex(sizeType globalIndex, int& nodeNum, sizeType& localIndex);

  /*
   * Description: Allocates a new ArrayChainNode to the ArrayChain 
   *
   * Return: True if allocation was successful
   *
   * Parameters: void
   */
  bool allocateNewNode();
};

#endif