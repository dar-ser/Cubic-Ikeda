#include "capd/capdlib.h"
#include <iomanip>
#include <iostream>
using namespace capd;
using namespace std;

void getAttractor(DTimeMap timeMap, DVector start,
                  double time, double skip, const string& filename);

void plotAttractor(const string &datafile, const string& outputfile);

