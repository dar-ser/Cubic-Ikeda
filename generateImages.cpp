#include "computations/pseudospectral.h"
#include "computations/newton.h"
#include "computations/pca.h"
#include "computations/cover.h"
#include "computations/pseudospectral.h"
#include "computations/bifurcation.h"
#include "computations/attractor.h"

#include <iomanip>
#include <iostream>
#include <algorithm>

using namespace capd;
using namespace std;
DMatrix approxMatrix;

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

int main()
{
    const int n = 7;
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

    /* PLOTTING APPROXIMATED ATTRACTOR */
    DVector x(n);
    for (int i = 1; i < n; i++)
        x[i] = 0.5;

    // making sure it's close to the attractor
    DTimeMap timeMap(solver);
    DTimeMap::SolutionCurve solution(0.);
    timeMap(200., x, solution);
    x = solution(200.);

    // plots approximated attractor
    double aStart = 1.5;
    double aChaotic = 1.542;
    double aEnd = 1.56;
    double aFrequency = 1000;
    int noSteps = 1000;

    CubicIkeda.setParameters({aStart});
    computeAttractorPoints(timeMap, x, 700., 200., "output/curve_start.csv");
    plotAttractor("output/curve_start.csv", "images/attractor_start.pdf");

    CubicIkeda.setParameters({aChaotic});
    computeAttractorPoints(timeMap, x, 700., 200., "output/curve_chaotic.csv");
    plotAttractor("output/curve_chaotic.csv", "images/attractor_chaotic.pdf");

    CubicIkeda.setParameters({aEnd});
    computeAttractorPoints(timeMap, x, 700., 200., "output/curve_end.csv");
    plotAttractor("output/curve_end.csv", "images/attractor_end.pdf");

    CubicIkeda.setParameters({a});

    /* PLOTTING BIFURCATION DIAGRAM */
    // warning: it will take a while
    computeBifurcationDiagram(CubicIkeda, pm, x, aStart, aEnd, aFrequency, noSteps, "output/bifurcationData.dat");
    plotBifurcationDiagram("output/bifurcationData.dat", "images/bifurcation.png", aStart, aEnd, aChaotic);


    /* SETTING UP FOR PCA AND ATTRACTOR IN NEW COORDINATES */
    DVector periodicPoint = getCandidate(CubicIkeda, pm, n, a, chaoticA);
    CubicIkeda.setParameters({chaoticA});
    
    DMatrix Dphi(n, n);
    DVector P = pm(periodicPoint, Dphi);
    DMatrix DP = pm.computeDP(P, Dphi);

    DMatrix DP_cut(n - 1, n - 1);
    for (int i = 0; i < n - 1; i++)
        for (int j = 0; j < n - 1; j++)
            DP_cut[i][j] = DP[i + 1][j + 1];

    DVector eigenRealPart(n - 1), eigenImPart(n - 1);
    DMatrix vectorRealPart(n - 1, n - 1), vectorImPart(n - 1, n - 1);
    alglib::computeEigenvaluesAndEigenvectors(DP_cut, eigenRealPart, eigenImPart,
                                              vectorRealPart, vectorImPart);

    std::map<double, DVector> eigenVectors;
    for (int i = 0; i < n - 1; i++)
        eigenVectors[eigenRealPart[i]] = vectorRealPart.column(i);

    std::sort(eigenRealPart.begin(), eigenRealPart.end(),
              [](double i, double j)
              { return abs(i) > abs(j); });

    interval data[] = {interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1)};
    IVector r(n, data);
    r *= 1e-5;
    r[0] *= 0;
    r[1] = interval(-1, 1) * 1e-6;

    IMatrix C(n, n);

    C[0][0] = 0;
    for (int i = 1; i < n; i++)
    {
        C[i][0] = 0;
        C[0][i] = 0;
        auto vec = eigenVectors[eigenRealPart[i - 1]];
        vec.normalize();
        for (int j = 1; j < n; j++)
            C[j][i] = vec[j - 1];
    }

    // plotting attractor in coordinates given by pca
    auto sample = samplePoints(C * r, periodicPoint);
    auto cloud = iteratePoincare(sample, pm, 100);
    auto V = mainDirections(cloud); // uses pca algorithm

    plotAttractor(pm, cloud, V, "images/attractor_pca.pdf");
}
