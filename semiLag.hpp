/* 
 * This header file contains basic
 * calculation function for dynamic
 * value advected by a incompressible
 * floid
 */

#ifndef _SEMI_LAGRANGE_
#define _SEMI_LAGRANGE_

#include <cmath>
#include <cstring>
#include <iostream>
#include "vector.hpp"

/* comment this out when
 * using it officially */
#include "myGrid.hpp"

#define GRID_SIZE particle_Grid_Len
/* This controlled deltaT between to grid */
#define TIME_STEP 0.10
/* This controlled the iteration in the semi-lagrangian scheme */
#define ITER_TIMES 40
/* This controlled the iterating time of the poisson solver */
#define POISSON_ITER 40

#define FETCH(i, j, k) ((i) * GRID_SIZE * GRID_SIZE + (j) * GRID_SIZE + (k))

/*
 *#define TEMPERATURE 0
 *#define DENSITY 1
 */

void semiLagrangeCalc(Grid &old, Grid &gen);

#endif
