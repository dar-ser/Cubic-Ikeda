#include "computations/pseudospectral.h"
#include "computations/bifurcation.h"
#include "computations/attractor.h"
#include <iomanip>
#include <iostream>
#include <algorithm>
using namespace capd;
using namespace std;

// the matrix needed for pseudospectral approximation
DMatrix approxMatrix;

/*
Pseudospectral approximation of Cubic Ikeda
f(x(t), x(t-tau)) = a(x(t-tau) - x(t-tau)^3)
*/

void approx_CubicIkeda(capd ::autodiff ::Node /*t*/, // unused time variable
                       capd ::autodiff ::Node in[],
                       int dimIn, // input variables x1 ,... , xn
                       capd ::autodiff ::Node out[],
                       int dimOut, // output : function values
                       capd ::autodiff ::Node params[],
                       int noParam // parameters
)
{
  int N = dimIn - 1;

  // out[0] = f(x(t), x(t-tau))
  //        = f(in[0], in[N])
  //        = a(in[N]-in[N]^3)
  capd ::autodiff ::Node a = params[0];
  capd ::autodiff ::Node xDelayed = in[N];
  out[0] = a * xDelayed * (1 - xDelayed * xDelayed);

  for (int i = 1; i <= N; i++)
  {
    capd ::autodiff ::Node result(0.);
    for (int j = 0; j <= N; j++)
      result += approxMatrix[i][j] * in[j];
    out[i] = result;
  }
}


/* MAIN FUNCTION */

int main(){
  const int N = 6;
  int dimIn = N + 1, dimOut = N + 1, noParams = 1, highestDerivative = 1;
  double a = 1.53;
  int taylorOrder = 20;
  double aChaotic = 1.542;
  double tau = 1;
  
  // for bifurcation diagram
  double aStart = 1.5;
  double aEnd = 1.56;
  double aFrequency = 1000;
  int noSteps = 1000;

  int n = N+1;
  string filenameM = "output/matrixM.txt";

  approxMatrix = compute_approxMatrix(tau, n, filenameM);
  cout << approxMatrix << endl << endl;

  DMap CubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
                  highestDerivative);
  CubicIkeda.setParameters({a});

  DOdeSolver solver(CubicIkeda, taylorOrder);
  DTimeMap timeMap(solver);

  // defining a Poincare map
  DCoordinateSection section(n, 0, 0);
  DPoincareMap pm(solver, section, poincare::MinusPlus);

  DVector x(n);
  for (int i = 1; i < n; i++)
    x[i] = 0.5;

  // making sure it's close to the attractor
  DTimeMap::SolutionCurve solution(0.);
  timeMap(200., x, solution);
  x = solution(200.); 

  // plots approximated attractor
  CubicIkeda.setParameters({aStart});
  getAttractor(timeMap, x, 700., 200., "output/curve_start.csv");
  plotAttractor( "output/curve_start.csv", "images/attractor_start.pdf");

  // CubicIkeda.setParameters({aChaotic});
  // getAttractor(timeMap, x, 700., 200., "output/curve_chaotic.csv");
  plotAttractor( "output/curve_chaotic.csv", "images/attractor_chaotic.pdf");

  // CubicIkeda.setParameters({aEnd});
  // getAttractor(timeMap, x, 700., 200., "output/curve_end.csv");
  plotAttractor( "output/curve_end.csv", "images/attractor_end.pdf");

  // CubicIkeda.setParameters({a});

  // plots bifurcation diagram: it will take a while
  // computeBifurcationDiagram(CubicIkeda, pm, x, aStart, aEnd, aFrequency, noSteps, "output/bifurcationData.dat");
  // plotBifurcationDiagram("output/bifurcationData.dat", "images/bifurcation.png", aStart, aEnd, aChaotic);
}