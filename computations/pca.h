#include <capd/capdlib.h>
#include <Eigen/Dense>
#include <vector>
#include <array>
#include <utility>
#include <numeric>
#include <iomanip>

using namespace capd;
using namespace Eigen;

DMatrix mainDirections(const std::vector<DVector> &pts);

std::vector<DVector> samplePoints(const IVector &rect, const DVector& x0);

std::vector<DVector> iteratePoincare(const std::vector<DVector> &sample,
                                    DPoincareMap &pm, int iterations);

void plotSection3D(const std::vector<DVector> &pts,
                   const DMatrix& V);

void plotColouredSection3D(DPoincareMap &pm, const std::vector<DVector> &pts, DMatrix& V);
// void plotDiagram(DPoincareMap &pm, const std::vector<DVector> &pts, DMatrix &V);
void plotDiagram(DPoincareMap &pm, 
                const std::vector<DVector> &pts,
                const DMatrix &V, 
                std::vector<IVector>& candidates);


void plotForLatex(DPoincareMap &pm,
                  const std::vector<DVector> &pts,
                  const DMatrix &V,
                  const std::string &filename);