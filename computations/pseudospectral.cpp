#include "pseudospectral.h"

/* functions needed for calculating approxMatrix */

DVector calculate_nodes(double tau, // time delay parameter
                        int n       // dimention
) {
  DVector s(n);
  const double step = M_PI / double(n-1);
  for (int i = 0; i < n; i++)
    s[i] = tau / 2. * (cos(double(i) * step) - 1.);
  return s;
}

DVector calculate_c(const DVector &s, // chebyshev nodes
                    int n             // dimention
) {
  DVector c(n);
  for (int i = 0; i < n; i++) {
    c[i] = 1.;
    for (int j = 0; j < n; j++) {
      if (j != i)
        c[i] *= (s[i] - s[j]);
    }
  }
  return c;
}

/* computing approxMatrix */

DMatrix compute_approxMatrix( double tau, int n, string filename) 
{
  DMatrix M(n, n);

  DVector s = calculate_nodes(tau, n);
  DVector c = calculate_c(s, n);

  for (int i = 0; i < n; i++) {
    double tmp = 0.;
    for (int j = 0; j < n; j++) {
      if (i != j) {
        double val = c[i] / (c[j] * (s[i] - s[j]));
        M[i][j] = val;
        tmp -= val;
      }
    }
    M[i][i] = tmp;
  }
  if (filename.compare("")) {
    std::ofstream out(filename);
    for (int i = 0; i < n; i++)
      out << M[i] << "\n";
    out.close();
  }
  return M;
}




void plot_chebyshev_nodes(int N, double tau,
                          const char* output = "images/chebyshev_nodes.pdf")
{
    FILE* gp = popen("gnuplot -persistent", "w");
    if (!gp)
        throw std::runtime_error("popen failed");

    // terminal & output
    fprintf(gp, "set terminal pdfcairo enhanced color font 'Palatino,15' size 14cm,9cm\n");
    fprintf(gp, "set output '%s'\n", output);

    // appearance
    fprintf(gp, "set style line 1 lc rgb '#2166AC' lw 2.2 dt 1\n");
    fprintf(gp, "set style line 2 lc rgb '#D73027' pt 7 ps 1.4\n");
    fprintf(gp, "set style line 3 lc rgb '#DDDDDD' lt 1 lw 0.4\n");

    fprintf(gp, "set grid ls 3\n");
    fprintf(gp, "set border lw 1.2\n");
    fprintf(gp, "set tics nomirror out\n");
    fprintf(gp, "set key off\n");

    fprintf(gp, "set xlabel 'k' font 'Palatino,15' offset 0,-0.5\n");
    fprintf(gp, "set ylabel 's_k' font 'Palatino,15' offset -0.5,0\n");

    // variables
    fprintf(gp, "tau = %.17g\n", tau);
    fprintf(gp, "N   = %d\n", N);
    fprintf(gp, "s(x) = tau/2.0 * (cos(pi*x/N) - 1.0)\n");
    fprintf(gp, "set xrange [-0.5 : N+0.5]\n");

    // data block
    std::string data;
    char buf[128];

    for (int k = 0; k <= N; ++k)
    {
        double sk = tau / 2.0 * (std::cos(k * M_PI / N) - 1.0);
        std::snprintf(buf, sizeof(buf), "%d %.17g\n", k, sk);
        data += buf;
    }

    fprintf(gp, "$nodes << EOD\n%sEOD\n", data.c_str());

    // labels
    fprintf(gp, "set label 1 's_0 = 0' at 0, s(0) left offset 1,-1 font 'Palatino,13' tc ls 2\n");
    fprintf(gp, "set label 2 '-{/Symbol t} = s_N' at N, s(N) right offset -0.4,1.25 font 'Palatino,13' tc ls 2\n");

    int id = 3;
    for (int k = 1; k < N; ++k)
    {
        double sk = tau / 2.0 * (std::cos(k * M_PI / N) - 1.0);

        double xoff = (k == 1 || k == N - 1) ? 0.0 : 0.5;
        double yoff = (k == 1) ? -1.6 : 1.3;

        fprintf(gp,
                "set label %d '%.4f' at %d, %.17g center offset %.1f, %.1f "
                "font 'Palatino,12' tc rgb '#555555'\n",
                id++, sk, k, sk, xoff, yoff);
    }

    // plot
    fprintf(gp,
            "plot [0:N] "
            "s(x) w l ls 1, "
            "$nodes u 1:2 w p ls 2\n");

    fflush(gp);
    pclose(gp);
}