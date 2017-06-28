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
/*
 *#include "myGrid.hpp"
 */

#define GRID_SIZE particle_Grid_Len
#define TIME_STEP 0.05
#define ITER_TIMES 10

/*
 *#define TEMPERATURE 0
 *#define DENSITY 1
 */

void semiLagrangeCalc(const Grid &old, Grid &gen);

#endif
