/*
 * File:  FieldVariable.h
 *
 * File Contents: Contains declarations for FieldVariable class
 *
 * Author: Yvonna
 * Contributors: Caroline Shung
 *               Nuttiiya Seekhao
 *               Kimberley Trickey
 */

#ifndef FIELDVARIABLE_H
#define FIELDVARIABLE_H

#include "../enums.h"

class World;

/*
 * FIELDVARIABLE CLASS DESCRIPTION:          Base class for per-patch field
 * data.
 */
class FieldVariable {
public:
  /*
   * Description:	Default FieldVariable constructor.
   *
   * Return: void
   *
   * Parameters: void
   */
  FieldVariable();

  /*
   * Description:	FieldVariable constructor.
   *
   * Return: void
   *
   * Parameters: orig  -- Pointer to an original FieldVariable
   */
  FieldVariable(const FieldVariable &orig);

  /*
   * Description:	Virtual FieldVariable destructor.
   *
   * Return: void
   *
   * Parameters: void
   */
  virtual ~FieldVariable();
};

#endif /* FIELDVARIABLE_H */