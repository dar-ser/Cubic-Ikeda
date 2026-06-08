#include "capd/capdlib.h"
#include <iomanip>
#include <iostream>
#include <cstdio>
#include <cmath>
#include <string>
#include <stdexcept>
using namespace capd;
using namespace std;

/* functions needed for calculating approxMatrix */

DVector calculate_nodes(double tau, // time delay parameter
                        int N       // dimention
);

DVector calculate_c(const DVector &s, // chebyshev nodes
                    int N             // dimention
);

/* computing approxMatrix */
/// @brief Computes matrix needed for approximating DDE with N-dimensional ODE
/// @param tau          time delay parameter
/// @param N            approximating vector's dimension - 1
/// @param filename     path to the output txt file
DMatrix compute_approxMatrix(double tau, int N, string filename);


void plot_chebyshev_nodes(int N, double tau, const char* output);