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
  // std::cout << " all: " << all << std::endl;
  return all;
}


void plotDiagram(DPoincareMap &pm, 
                const std::vector<DVector> &pts,
                const DMatrix &V, 
                std::vector<IVector>& candidates) {

  std::cout << " plotting 'diagram' " << std::endl;

  const int n = pts.size();
  const int dim = pts[0].dimension();

  DVector mean(dim);
  for (auto &p : pts)
    mean += p;
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

  double minLeftVal = 1e99;
  double minRightVal = 1e99;
  double maxLeftVal = 1e99;
  double maxRightVal = 1e99;
  int indexMinLeft = -1;
  int indexMinRight = -1;
  int indexMaxLeft = -1;
  int indexMaxRight = -1;
  double leftMin = -0.0161;
  double leftMax = 0.03;
  double rightMin = 0.035;
  double rightMax = 0.099;

  std::vector<DVector> id;

  std::vector<DVector> tp(n), ti(n);
  for (int i = 0; i < n; ++i) {
    tp[i] = project(pts[order[i]]);
    ti[i] = project(pm(pm(pts[order[i]])));

    if(std::abs(tp[i][0] - ti[i][0]) < 1e-4){
      id.push_back(tp[i]);
    }


    double tmp = std::abs(tp[i][0] - leftMin);
    if (tmp < minLeftVal) {
      indexMinLeft = i;
      minLeftVal = tmp;
      continue;
    }
    tmp = std::abs(tp[i][0] - leftMax);
    if (tmp < maxLeftVal) {
      indexMaxLeft = i;
      maxLeftVal = tmp;
      continue;
    }
    tmp = std::abs(tp[i][0] - rightMin);
    if (tmp < minRightVal) {
      indexMinRight = i;
      minRightVal = tmp;
      continue;
    }
    tmp = std::abs(tp[i][0] - rightMax);
    if (tmp < maxRightVal) {
      indexMaxRight = i;
      maxRightVal = tmp;
      continue;
    }
  }


  // candidates to check covering relation
  IVector left(dim), right(dim);
  for(int i = 0; i < dim; i++){
    left[i] = interval(tp[indexMinLeft][i], tp[indexMaxLeft][i]);
    right[i] = interval(tp[indexMinRight][i], tp[indexMaxRight][i]);
  }

  candidates[0]=left;
  candidates[1]=right;

  /// finding a point near the end of the tail
  // to see where it goes
  std::vector<double> v = {-0.012,-0.005,0.0005};
  double min = 1e99;
  int index = -1;
  for(int i = 0; i < n; i++){
    double tmp = std::abs(tp[i][0] - v[0]) + std::abs(tp[i][1] - v[1]) + std::abs(tp[i][2] - v[2]);
    if (tmp < min) {
      index = i;
      min = tmp;
    }
  }


  // the largest values across all three axes
  double lo[2] = {2147483647, 2147483647};
  double hi[2] = {-2147483647, -2147483647};
  for (int i = 0; i < n; ++i) {
    auto p0 = tp[i];
    auto p1 = ti[i];
    lo[0] = std::min(lo[0], p0[0]);
    lo[1] = std::min(lo[1], p1[0]);
    hi[0] = std::max(hi[0], p0[0]);
    hi[1] = std::max(hi[1], p1[0]);
  }


  auto maxSpan = std::max(hi[0] - lo[0],hi[1] - lo[1]);
  double h = maxSpan / 8;

  std::cout << "here\n";
  FILE *gp = popen("gnuplot -persistent", "w");
  fprintf(gp, "set terminal qt enhanced\n");
  fprintf(gp, "set palette defined (0 'blue',1 'cyan',2 'green',3 'yellow',4 "
               "'red',5 'magenta')\n");
  fprintf(gp, "set cbrange [0:1]\n");
  fprintf(gp, "unset colorbox\n");
  fprintf(gp, "set mouse\n");
  fprintf(gp, "set border 4095\n");
  fprintf(gp, "set title 'P(PCA)'\n");
  fprintf(gp, "set xlabel 'e0'\nset ylabel 'e0'\n");
  fprintf(gp, "set xrange [%g:%g]\nset yrange [%g:%g]\n",
          lo[0] - h, hi[0] + h, lo[1] - h, hi[1] + h);
  fprintf(gp, "plot '-' u 1:2:3 w points pt 7 ps 0.3 lc palette notitle, \
              '-' w points pt 7 ps 2 lc rgb 'red' title 'fixed points', \
              '-' w points pt 7 ps 2 lc rgb 'black' title 'left set', \
              '-' w points pt 7 ps 2 lc rgb 'green' title 'right set'\n");
  for (int i = 0; i < n; ++i)
  fprintf(gp, "%g %g %g\n", tp[i][0], ti[i][0], i / (double)(n - 1));
  fprintf(gp, "e\n");

  std::cout << "left min : " << tp[indexMinLeft][0] <<  " -> " << ti[indexMinLeft][0] << std::endl;
  std::cout << "left min : " << tp[indexMaxLeft][0] <<  " -> " << ti[indexMaxLeft][0] << std::endl;
  std::cout << "right max : " << tp[indexMinRight][0] <<  " -> " << ti[indexMinRight][0] << std::endl;
  std::cout << "right max : " << tp[indexMaxRight][0] <<  " -> " << ti[indexMaxRight][0] << std::endl;
  
  for(auto& v : id){
    fprintf(gp, "%g %g\n", v[0], v[0]);
  }
  fprintf(gp, "e\n");

  fprintf(gp, "%g %g\n", tp[indexMinLeft][0], ti[indexMinLeft][0]);
  fprintf(gp, "%g %g\n", tp[indexMaxLeft][0], ti[indexMaxLeft][0]);
  fprintf(gp, "e\n");

  fprintf(gp, "%g %g\n", tp[indexMinRight][0], ti[indexMinRight][0]);
  fprintf(gp, "%g %g\n", tp[indexMaxRight][0], ti[indexMaxRight][0]);
  fprintf(gp, "e\n");
  
  fprintf(gp, "pause mouse close\n");
  pclose(gp);
}


void plotForLatex(DPoincareMap &pm,
                  const std::vector<DVector> &pts,
                  const DMatrix &V,
                  const std::string &filename)
{
    std::cout << "plotting attractor (latex pdf)" << std::endl;

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