// #include "computations/pseudospectral.h"
// #include "computations/newton.h"
// #include "computations/pca.h"
// #include "computations/cover.h"

// #include <iomanip>
// #include <iostream>
// #include <algorithm>

// using namespace capd;
// using namespace std;

// DMatrix approxMatrix;

// void approx_CubicIkeda(capd ::autodiff ::Node /*t*/, // unused time variable
//                     capd ::autodiff ::Node in[],
//                     int dimIn, // input variables x1 ,... , xn
//                     capd ::autodiff ::Node out[],
//                     int dimOut, // output : function values
//                     capd ::autodiff ::Node params[],
//                     int noParam // parameters
// )
// {
//     int N = dimIn - 1;

//     capd ::autodiff ::Node a = params[0];
//     capd ::autodiff ::Node xDelayed = in[N];
//     out[0] = a * xDelayed * (1 - xDelayed * xDelayed);

//     for (int i = 1; i <= N; i++)
//     {
//         capd ::autodiff ::Node result(0.);
//         for (int j = 0; j <= N; j++)
//             result += approxMatrix[i][j] * in[j];
//         out[i] = result;
//     }
// }


// int main()
// {
//     int n = 7;
//     int dimIn = n, dimOut = n, noParams = 1, highestDerivative = 1;
//     double a = 1.53;
//     int taylorOrder = 20;
//     double returnTime = 0.;
//     double chaoticA = 1.542;
//     double tau = 1;

//     string filename = "";

//     approxMatrix = compute_approxMatrix(tau, n, filename);

//     DMap CubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
//                     highestDerivative);
//     CubicIkeda.setParameters({a});

//     DOdeSolver solver(CubicIkeda, taylorOrder);

//     /* defining Poincare map */
//     DCoordinateSection section(n, 0, 0);
//     DPoincareMap pm(solver, section, poincare::MinusPlus);
//     DVector periodicPoint = getCandidate(CubicIkeda, pm, n, a, chaoticA);
//     DMatrix Dphi(n, n); 
//     DVector P = pm(periodicPoint, Dphi);  
//     DMatrix DP = pm.computeDP(P, Dphi);
//     DMatrix DP_cut(n-1, n-1);
//     for (int i = 0; i < n-1; i++)
//         for (int j = 0; j < n-1; j++)
//             DP_cut[i][j] = DP[i + 1][j + 1];

//     DVector eigenRealPart(n-1), eigenImPart(n-1);
//     DMatrix vectorRealPart(n-1, n-1), vectorImPart(n-1, n-1);
//     alglib::computeEigenvaluesAndEigenvectors(DP_cut, eigenRealPart, eigenImPart,
//                                               vectorRealPart, vectorImPart);

//     std::map<double, DVector> eigenVectors;
//     for (int i = 0; i < n-1; i++)
//         eigenVectors[eigenRealPart[i]] = vectorRealPart.column(i);

//     std::sort(eigenRealPart.begin(), eigenRealPart.end(),
//               [](double i, double j){ return abs(i) > abs(j); });

//     interval data[] = {interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1)};
//     IVector r(n, data);

//     r *= 1e-5;
//     r[0] *= 0;
//     r[1] = interval(-1, 1) * 1e-6;

//     IMap ICubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
//                      highestDerivative);
//     ICubicIkeda.setParameters({chaoticA});

//     IVector x0(periodicPoint);
//     x0[0] = 0.;
//     IMatrix B(n, n);
//     auto fx0 = ICubicIkeda(x0);
//     fx0.normalize();

//     B[0][0] = fx0[0];
//     for (int i = 1; i < n; i++)
//     {
//         B[i][0] = fx0[i];
//         B[0][i] = 0;
//         auto vec = eigenVectors[eigenRealPart[i - 1]];
//         vec.normalize();
//         for (int j = 1; j < n; j++)
//             B[j][i] = vec[j-1];
//     }

//     IMatrix A = matrixAlgorithms::gaussInverseMatrix(B);

//     for (int i = 0; i < n; i++)
//         B[i][0] = 0;

//     auto sample  = samplePoints(B * r, periodicPoint);
//     auto cloud = iteratePoincare(sample, pm, 100);

//     auto V = mainDirections(cloud);
//     plotColouredSection3D(pm, cloud, V);

//     // IVector first(n), second(n);
//     // plotForLatex(pm,cloud,V, "images/attractor_pca.pdf");

// }
#include "computations/pseudospectral.h"
#include "computations/newton.h"
#include "computations/pca.h"
#include "computations/cover.h"
#include "capd/covrel/HSetMD.hpp"
#include "capd/dynsys/DynSysMap.h"

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
    DVector periodicPoint = getCandidate(CubicIkeda, pm, n, a, chaoticA);
    CubicIkeda.setParameters({chaoticA}); 
    DMatrix Dphi(n, n);
    DVector P = pm(periodicPoint, Dphi);
    DMatrix DP = pm.computeDP(P, Dphi);

    DMatrix DP_cut(n-1,n-1);
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

    interval data[] = {interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1), interval(-1, 1)};
    IVector r(n, data);
    r *= 1e-5;
    r[0] *= 0;
    r[1] = interval(-1, 1) * 1e-6;


    IMap ICubicIkeda(approx_CubicIkeda, dimIn, dimOut, noParams,
                    highestDerivative);
    ICubicIkeda.setParameters({chaoticA});

    IVector x0(periodicPoint);
    x0[0] = 0.;

    IMatrix C(n,n);
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
    for (int i = 0; i < n; i++)
        C[i][0] = 0;

    

    auto sample  = samplePoints(C * r, periodicPoint);
    auto cloud = iteratePoincare(sample, pm, 100);


    auto V = mainDirections(cloud);

    plotForLatex(pm, cloud, V, "images/attractor_pca.pdf");


}


