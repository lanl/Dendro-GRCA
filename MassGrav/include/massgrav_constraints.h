#ifndef MASSGRAV_CONSTRAINTS_H
#define MASSGRAV_CONSTRAINTS_H

#include <iostream>
#include "parameters.h"
#include "grUtils.h"

void enforce_massgrav_constraints(double **uiVar, const unsigned int offset);

#endif
