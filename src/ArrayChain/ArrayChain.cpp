/*
 * ArrayChain.h
 *
 * File Contents: Contains the ArrayChain class.
 *
 * Author: NungnunG
 * Contributors: Caroline Shung
 *               Kimberley Trickey
 */

#include "ArrayChain.h"

template<typename T>
ArrayChain<T>::ArrayChain() {
	head = NULL;
	tail = NULL;
	numNodes = 0;
	chainCapacity = 0;
	nodeCapacity = 0;
}

template<typename T>
ArrayChain<T>::ArrayChain(sizeType nodeCap, sizeType chainCap, PrintElementFunc pFunc, DestroyElementFunc dFunc) {
	head = NULL;
	tail = NULL;
	numNodes = 0;
	chainCapacity = chainCap;
	nodeCapacity = nodeCap;
	printFunc = pFunc;
	destroyFunc = dFunc;
}

template<typename T>
ArrayChain<T>::ArrayChain(const ArrayChain<T>& source) {
	head = source.head;
	tail = source.tail;
	numNodes = source.numNodes;
	chainCapacity = source.chainCapacity;
	nodeCapacity = source.nodeCapacity;
	printFunc = source.printFunc;
	destroyFunc = source.destroyFunc;
}

template<typename T>
ArrayChain<T>& ArrayChain<T>::operator=(const ArrayChain<T>& source) {
	if (this == &source) return *this;

	head = source.head;
	tail = source.tail;
	numNodes = source.numNodes;
	chainCapacity = source.chainCapacity;
	nodeCapacity = source.nodeCapacity;
	printFunc = source.printFunc;
	destroyFunc = source.destroyFunc;
	return *this;
}

template<typename T>
ArrayChain<T>::~ArrayChain() {
	// Traverse through the chain and destroy all things stored in nodes
	ArrayChainNode<T> *currNode = head;
	ArrayChainNode<T> *prevNode = NULL;
	while (currNode != NULL) {
		prevNode = currNode;
		currNode = currNode->next;
		delete prevNode;
	}
}

template<typename T>
bool ArrayChain<T>::addData(T val, int tid) {
	// Traverse through the chain to find available element in a node
	ArrayChainNode<T> *currNode = head;
	while (currNode != NULL) {
		if (currNode->addData(val, tid)) return true;
		currNode = currNode->next;
	}
	// If all nodes are full, allocate a new node and add data there
	if (!allocateNewNode()) return false;
	return tail->addData(val, tid);
}

template<typename T>
bool ArrayChain<T>::deleteData(sizeType globalIndex, int tid) {
	int nodeNum;
	sizeType localIndex;
	ArrayChainNode<T> *currNode;
	//  std::cerr << "	Deleting cell " << globalIndex << " from chain" << std::endl;
	if (!head) {  // Empty ArrayChain
		std::cout << "Warning: deleteData() trying to delete from empty chain" << std::endl;
		return false;
	}
	// Convert global to local index
	globalToLocalIndex(globalIndex, nodeNum, localIndex);
	// Traverse the chain to get to nodeNum-th node and delete the data there
	currNode = head;
	for (int node = 0; node < nodeNum; node++) currNode = currNode->next;
	return currNode->deleteData(localIndex, tid);
}

template<typename T>
bool ArrayChain<T>::concAddData(T val, int tid) {
	std::cout << "ArrayChain<T>::concAddData() not implemented yet" << std::endl;
	return false;
}

template<typename T>
T& ArrayChain<T>::getDataAt(sizeType globalIndex) {
	int nodeNum;
	sizeType localIndex;
	ArrayChainNode<T> *currNode;
	// Convert global to local index
	globalToLocalIndex(globalIndex, nodeNum, localIndex);
	// Traverse the chain to get to nodeNum-th node and get the data there
	currNode = head;
	for (int node = 0; node < nodeNum; node++) currNode = currNode->next;
	return currNode->getDataAt(localIndex);
}

template<typename T>
sizeType ArrayChain<T>::size() {
	if (!numNodes) return 0;
	sizeType sz = (numNodes * nodeCapacity) - tail->getNumFreeApprox();
	return sz;
}

template<typename T>
sizeType ArrayChain<T>::actualSize() {
	sizeType sz = 0;
	ArrayChainNode<T> *curr = head;
	while (curr) {
		sz += curr->actualSize();
		curr = curr->next;
	}
	return sz;
}

template<typename T>
void ArrayChain<T>::print() {
	sizeType numElementsTotal = numNodes * nodeCapacity;
	std::cout << "[";
	if (printFunc) {
		int i = 0;
		for (i = 0; i < (numElementsTotal - 1); i++) {
			std::cout << " ";
			printFunc(getDataAt(i));
			std::cout << ",";
		}
		std::cout << " ";
		printFunc(getDataAt(i));
		std::cout << "]" << std::endl;
	} else {
		std::cout << "ArrayChain<T>::print(): Print function not yet registered!" << std::endl;
		//    int i = 0;
		//    for (i = 0; i < (numElementsTotal - 1); i++) {
		//      std::cout << " " << getDataAt(i) << ",";
		//    }
		//    std::cout << " " << getDataAt(i) << "]" << std::endl;
	}
}

template<typename T>
void ArrayChain<T>::globalToLocalIndex(sizeType globalIndex, int& nodeNum, sizeType& localIndex) {
	nodeNum = globalIndex / nodeCapacity;
	localIndex = globalIndex % nodeCapacity;
}

template<typename T>
bool ArrayChain<T>::allocateNewNode() {
	ArrayChainNode<T>* newNode = new ArrayChainNode<T>(nodeCapacity, destroyFunc);
	std::cout << "  allocating new node" << std::endl;
	std::cerr << "  allocating new node" << std::endl;
	if (!newNode) {
		return false;
	}
	if (!head) {  // Empty ArrayChain
		head = newNode;
		tail = head;
	} else {
		tail->next = newNode;
		tail = newNode;
	}
	numNodes++;
	return true;
}

template class ArrayChain<int *>;
template class ArrayChain<Cell *>;