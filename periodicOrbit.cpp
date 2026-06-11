#include "computations/pseudospectral.h"
#include "computations/newton.h"
#include "computations/cover.h"
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
                       int dimIn, // input variables x0 ,... , xN
                       capd ::autodiff ::Node out[],
                       int dimOut, // output : function values
                       capd ::autodiff ::Node params[],
                       int noParam // parameters
)
{
    int N = dimIn - 1;

    // out[0] = f(x(t), x(t-tau))
    //        = f(in[0], in[N])
    //        = a(in[N]-in[N])^3
    capd ::autodiff ::Node a = params[0];
    capd ::autodiff ::Node xDelayed = in[N];
    out[0] = a * xDelayed * (1 - xDelayed * xDelayed);

    for (int i = 1; i < dimIn; i++)
    {
        capd ::autodiff ::Node result(0.);
        for (int j = 0; j <= N; j++)
            result += approxMatrix[i][j] * in[j];
        out[i] = result;
    }
}

/* MAIN FUNCTION */

int main()
{
    int n = 7;
    int dimIn = n, dimOut = n, noParams = 1, highestDerivative = 1;
    double a = 1.53;
    int taylorOrder = 20;
    double returnTime = 0.;
    double chaoticA = 1.542;
    double tau = 1;

    string filename = "";

    approxMatrix = compute_approxMatrix(tau, n, filename);

    DMap CubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
                    highestDerivative);
    CubicIkeda.setParameters({a});

    DOdeSolver solver(CubicIkeda, taylorOrder);

    /* defining Poincare map */
    DCoordinateSection section(n, 0, 0);
    DPoincareMap pm(solver, section, poincare::MinusPlus);


    // defining a rough approximation of a periodic orbit
    DVector periodicPoint = getCandidate(CubicIkeda, pm, n, a, chaoticA);
    

    // setting up new coordinate system around that point for interval set definition

    // getting eigenvectors and eigenvalues of DP(periodicPoint) 
    // without the parts corresponding to the first coordinate
    // (the section enforces it to be zero)

    DMatrix Dphi(n, n); 
    DVector P = pm(periodicPoint, Dphi); 
    // Dphi contains the derivative of the flow evaluated at (periodicPoint,returnTime)
 
    DMatrix DP = pm.computeDP(P, Dphi);
    DMatrix DP_cut(n-1, n-1);
    for (int i = 0; i < n-1; i++)
        for (int j = 0; j < n-1; j++)
            DP_cut[i][j] = DP[i + 1][j + 1];

    DVector eigenRealPart(n-1), eigenImPart(n-1);
    DMatrix vectorRealPart(n-1, n-1), vectorImPart(n-1, n-1);
    alglib::computeEigenvaluesAndEigenvectors(DP_cut, eigenRealPart, eigenImPart,
                                              vectorRealPart, vectorImPart);

    std::map<double, DVector> eigenVectors;
    for (int i = 0; i < n-1; i++)
        eigenVectors[eigenRealPart[i]] = vectorRealPart.column(i);

    std::sort(eigenRealPart.begin(), eigenRealPart.end(),
              [](double i, double j){ return abs(i) > abs(j); });

    cout << endl
         << "real part of eigenvalues:\n"
         << eigenRealPart << endl
         << "imaginary part of eigenvalues:\n"
         << eigenImPart << endl
         << endl;
    /* so we see that we only have real eigenvalues :) */

    
    // definition of the interval set centred around periodicPoint
    // is given by  x0 + C * r
    // where C - the basis transformation matrix for the set
    //   and A - the basis transformation matrix for the Poincare map's output

    // we define r specifically for n = 7 
    interval data[] = {interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1)};
    IVector r(n, data);

    r *= 1e-6;
    r[0] *= 0;
    r[1] = interval(-1, 1) * 1e-7;

    IMap ICubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
                     highestDerivative);
    ICubicIkeda.setParameters({chaoticA});

    IVector x0(periodicPoint);
    x0[0] = 0.;

    // the basis transformation matrices are defined using the eigenvectors computed before 
    // with f(x0)/||f(x0)|| in the first column for A^{-1}
    // and zeros in the first column for C
    IMatrix C(n, n);
    auto fx0 = ICubicIkeda(x0);
    fx0.normalize();

    C[0][0] = fx0[0];
    for (int i = 1; i < n; i++)
    {
        C[i][0] = fx0[i];
        C[0][i] = 0;
        auto vec = eigenVectors[eigenRealPart[i - 1]];
        vec.normalize();
        for (int j = 1; j < n; j++)
            C[j][i] = vec[j-1];
    }

    IMatrix A = matrixAlgorithms::gaussInverseMatrix(C);

    // for more precise computation we ensure that the set is contained in the section
    for (int i = 0; i < n; i++)
        C[i][0] = 0;

    cout << "B:\n" << C << endl
         << "A:\n" << A << endl;

    // defining the set s
    C1Rect2Set s(x0, C, r);

    // defining interval Poincare map
    IOdeSolver ISolver(ICubicIkeda, taylorOrder);
    ICoordinateSection ISection(n, 0, 0);
     
    IPoincareMap IPm(ISolver, ISection, poincare::MinusPlus);
    interval IReturnTime = 0;

    // computing P(s) in coordinates given by A
    IVector P0 = IPm(s, x0, A, IReturnTime);
    
    // we know that P0 contains the set P(r) 
    cout << endl << "r: " << r << endl;
    cout << "P0: " << P0 << endl << endl;

    cout << "checking covering relations:\n";
    cout << std::boolalpha;

    bool cover = checkSelfCover(IPm, A, C, r, x0);

    cout << "covering : " << cover << endl;
}
