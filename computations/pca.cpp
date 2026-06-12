#include "pca.h"

// "Sorts" directions of a DVector point cloud
// Returns matrix that contains new vectors spanning the space 
//    sorted by descending variance.
DMatrix mainDirections(const std::vector<DVector> &pts) {
  std::cout << " finding main directions \n";

  const int n = (int)pts.size();
  const int dim = (int)pts[0].dimension();

  MatrixXd X(n, dim);
  for (int i = 0; i < n; ++i)
    for (int k = 0; k < dim; ++k)
      X(i, k) = pts[i][k];

  // computing means of every coordinate,
  // then subtracting them from each point
  VectorXd mean = X.colwise().mean();
  X.rowwise() -= mean.transpose();

  MatrixXd C = (X.transpose() * X) / (n - 1); // covariance

  SelfAdjointEigenSolver<MatrixXd> eig(C);
  MatrixXd V = eig.eigenvectors().rowwise().reverse();
  // this way V.col(0) is the direction of maximum value

  // transform it to DMatrix
  DMatrix V1(dim, dim);
  for (int i = 0; i < dim; i++)
    for (int j = 0; j < dim; j++)
      V1[i][j] = V(i, j);
  V1 = matrixAlgorithms::gaussInverseMatrix(V1);
  return V1;
}

std::vector<DVector> samplePoints(const IVector &rect, const DVector &x0) {
  std::cout << " generating sample points \n";

  std::vector<DVector> pts;
  std::vector<double> coord1;

  double low = rect[1].leftBound();
  double mid = rect[1].mid().leftBound();
  double length = mid - low;

  while (low < mid) {
    coord1.push_back(low);
    low += length * .01;
  }

  for (auto first : coord1) {
    DVector tmp = x0;
    tmp[1] = first;
    pts.push_back(tmp);
  }

  return pts;
}

std::vector<DVector> iteratePoincare(const std::vector<DVector> &pts,
                                     DPoincareMap &pm, int iterations) {
  std::cout << " iterating sample points \n";

  std::vector<DVector> all;
  all.reserve(pts.size() * iterations);
  for (size_t i = 0; i < pts.size(); ++i) {
    DVector p = pts[i];
    for (int j = 0; j < iterations; ++j) {
      p = pm(p);
      all.push_back(p);
    }
  }
  return all;
}

void plotAttractor(DPoincareMap &pm,
                  const std::vector<DVector> &pts,
                  const DMatrix &V,
                  const std::string &filename)
{
    std::cout << " plotting attractor" << std::endl;

    const int n = pts.size();
    const int dim = pts[0].dimension();

    DVector mean(dim);
    for (auto &p : pts) mean += p;
    mean /= (double)n;

    auto project = [&](const DVector &p) {
        DVector c(dim);
        for (int k = 0; k < dim; ++k)
            c[k] = p[k] - mean[k];
        return V * c;
    };

    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b) {
        return project(pts[a])[0] < project(pts[b])[0];
    });

    std::vector<DVector> tp(n), ti(n);
    for (int i = 0; i < n; ++i) {
        tp[i] = project(pts[order[i]]);
        ti[i] = project(pm(pts[order[i]]));
    }

    double lo[3] = { 2147483647,  2147483647,  2147483647};
    double hi[3] = {-2147483647, -2147483647, -2147483647};
    for (int i = 0; i < n; ++i)
        for (int a = 0; a < 3; ++a) {
            lo[a] = std::min(lo[a], ti[i][a]);
            hi[a] = std::max(hi[a], ti[i][a]);
        }

    double maxSpan = 0;
    for (int a = 0; a < 3; ++a)
        maxSpan = std::max(maxSpan, hi[a] - lo[a]);

    double cx = (lo[0]+hi[0])/2, cy = (lo[1]+hi[1])/2, cz = (lo[2]+hi[2])/2;
    double h = maxSpan / 2;

    const double rot_x = 60.0;
    const double rot_z = 30.0;

    FILE *gp = popen("gnuplot", "w");

    fprintf(gp, "set terminal pdfcairo enhanced color size 6,6\n");
    fprintf(gp, "set output '%s'\n", filename.c_str());

    fprintf(gp, "set palette defined (0 'blue',1 'cyan',2 'green',3 'yellow',4 'red',5 'magenta')\n");
    fprintf(gp, "set cbrange [0:1]\n");
    fprintf(gp, "unset colorbox\n");
    fprintf(gp, "unset key\n");
    fprintf(gp, "set border 4095\n");
    fprintf(gp, "set border lc rgb '#808080'\n");

    fprintf(gp, "set tics font ',14'\n");
    fprintf(gp, "set tics out\n");
    fprintf(gp, "set tics scale 1.0\n");
    fprintf(gp, "set xtics auto offset 0,-1\n");
    fprintf(gp, "set ytics auto offset 0,-1\n");
    fprintf(gp, "set ztics auto\n");

    fprintf(gp, "set xlabel 'e0' font ',22' offset 0,-2\n");
    fprintf(gp, "set ylabel 'e1' font ',22' offset 7.5,0\n");
    fprintf(gp, "set zlabel 'e2' font ',22' offset 0.75,0\n");

    fprintf(gp, "set view %g, %g, 1, 1\n", rot_x, rot_z);

    fprintf(gp, "set xrange [%g:%g]\n", cx-h, cx+h);
    fprintf(gp, "set yrange [%g:%g]\n", cy-h, cy+h);
    fprintf(gp, "set zrange [%g:%g]\n", cz-h, cz+h);
    fprintf(gp, "set xyplane at %g\n", cz-h);

    fprintf(gp, "splot '-' u 1:2:3:4 w points pt 7 ps 0.25 lc palette notitle\n");

    for (int i = 0; i < n; ++i)
        fprintf(gp, "%g %g %g %g\n",
                ti[i][0], ti[i][1], ti[i][2],
                i / (double)(n - 1));

    fprintf(gp, "e\n");
    fflush(gp);
    pclose(gp);
}