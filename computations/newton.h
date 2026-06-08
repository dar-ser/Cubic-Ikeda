#include "capd/capdlib.h"
#include <iomanip>
#include <iostream>
using namespace capd;
using namespace std;


/// @brief Uses Newton's method for finding fixed point of the function map for the parameter aEnd,
///        by atarting at the lower value aStart and moving values higher.
/// @param dmap     map used for integration.
/// @param pm       Poincaré map defined on the solver of dmap and a section.
/// @param n        dimension of the map's arguments
/// @param aStart   starting, lower value of the parameter a.
/// @param aEnd     ending, higher (so more chaotic) value of the parameter a.
DVector getCandidate(DMap& map, DPoincareMap& pm, int n, double aStart, double aEnd);
