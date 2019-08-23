// all.h is a header file wich includes all the other headers we need in this project. Both standard c/c++ libraries
// and also our libraries such as ALGLIB

#ifndef ALL_H
#define ALL_H

// General includes:

//#include <stdio.h>
//#include <math.h>
#include <iostream>
#include <fstream>
#include <time.h>
#include <omp.h>
#include <cmath>

// Our includes:

#include "point.h"
#include "globals.h"
#include "setup.h"
#include "misc_calc.h"
#include "iterative.h"
#include "iterative_exp.h"
#include "alglibinternal.h"
#include "alglibmisc.h"
#include "integration.h"

// Constants:

double const pi = 3.14159265;
double const pisqrt = 1.772453;
double const e = 2.71828183;

double const I_inc_1 = 0.0;
double const I_inc_2 = 0.0;
double const I_inc_3 = 1.0;
double const I_inc_4 = 1.0;

int const periodic_boundary = 1;


#endif

