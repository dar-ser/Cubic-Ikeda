#include "newton.h"

// A helper function implementing the Newton's method for finding a zero of the function pm
// starting at the point x0 with precision given as a parameter
DVector getZero(DPoincareMap &pm, DVector &x0, double precision)
{
    DVector xn(x0);
    int n = xn.dimension();

    DMatrix Dphi(n, n);
    DMatrix Id = vectalg::Matrix<double, 0, 0>::Identity(n);

    DVector P = pm(xn, Dphi); 
    // Dphi contains the derivative of the flow evaluated at (xn,returnTime(xn))

    // but we actually care about the function P-Id
    DMatrix DF = pm.computeDP(P, Dphi) - Id;
    DMatrix invDF = matrixAlgorithms::gaussInverseMatrix(DF);

    int i = 0;
    capd::vectalg::MaxNorm<DVector, DMatrix> maxNorm;
    while (maxNorm(P - xn) > precision)
    {
        i++;
        xn = xn - invDF * (P - xn);
        P = pm(xn, Dphi);
        DF = pm.computeDP(P, Dphi) - Id;
        invDF = matrixAlgorithms::gaussInverseMatrix(DF);
    }

    return xn;
}


DVector getCandidate(DMap &map, DPoincareMap &pm, int n, double aStart, double aEnd)
{
    // searching first for a fixed point for the lower (so less chaotic) parameter a = aStart
    auto entry_a = map.getParameter(0);
    map.setParameters({aStart});
    double precision = 1e-15;

    DVector start(n);
    start[0] = 0;
    start[n-1] = 1.44;
    for (int i = 1; i < n-1; i++)
        start[i] = 1 / 4.;

    DVector candidate = getZero(pm, start, precision);

    // the point is to look for fixed points higher (so more chaotic) parameter a values
    // by starting at the last calculated fixed point, which is close to the new fixed point 
    double a = aStart;
    while (a < aEnd)
    {
        a += .002;
        map.setParameters({a});
        candidate = getZero(pm, candidate, precision); 
    }

    // calculating the fixed point for the parameter value we are interested in
    map.setParameters({aEnd});
    candidate = getZero(pm, candidate, precision); 

    map.setParameters({entry_a});
    cout << candidate << " " << pm(candidate) << endl;
    return candidate;
}

