#include "capd/capdlib.h"
#include <iomanip>
#include <iostream>
using namespace capd;
using namespace std;

/// @brief Uses gnuplot to plot a bifurcation diagram for the Poincaré map
///        of the cubic Ikeda function.
/// @param dmap     map used for integration.
/// @param pm       Poincaré map defined on the solver of dmap and a section.
/// @param x        initial state vector.
/// @param aStart   starting value of parameter a.
/// @param aEnd     ending value of parameter a.
/// @param noSteps  number of iterations computed for each value of a.
/// @param chaoticA specific value to plot.  
void plotBifurcationDiagram(DMap &dmap, DPoincareMap &pm, const DVector &x,
                            double aStart, double aEnd, double aFrequency,
                            int noSteps, double chaoticA);



/// @brief Computes bifurcation diagram data and saves it to a file for later plotting.
/// @param CubicIkeda   map used for integration.
/// @param pm           Poincaré map defined on the solver of dmap and a section.
/// @param x            initial state vector.
/// @param aStart       start of parameter range.
/// @param aEnd         end of parameter range.
/// @param aFrequency   number of parameter sampling steps.
/// @param noSteps      number of iterations per parameter value.
/// @param filename     output file for (a, x) data.
void computeBifurcationDiagram(DMap &CubicIkeda, DPoincareMap &pm, const DVector &x,
                                double aStart, double aEnd, double aFrequency,
                                int noSteps, const string &filename);

/// @brief Plots a bifurcation diagram from a precomputed data file using Gnuplot.
/// @param datafile     input file containing (a, x) data.
/// @param outputfile   output image file (png/pdf)
/// @param x            initial state vector.
/// @param aStart       start of parameter range.
/// @param aEnd         end of parameter range.
/// @param chaoticA     specific parameter a value for highlighted in the plot.
void plotBifurcationDiagram(const string &datafile, const string &outputfile,
                            double aStart, double aEnd, double chaoticA);