#include "attractor.h"

void getAttractor(DTimeMap timeMap, DVector start,
                      double time, double skip, const string& filename) {
  DTimeMap::SolutionCurve solution(0.);
  timeMap(skip, start, solution);
  start = solution(skip); 
  
  timeMap(time, start, solution);

  int N = start.dimension() - 1;

  if (filename.compare("")) {
    ofstream csv(filename);
    csv << "x_present,x_delayed\n";
    for (double t = 0; t <= time - skip; t += 0.05) {
      DVector point = solution(t);
      csv << point[0] << "," << point[N] << "\n";
    }
    csv.close();
  }
}

void plotAttractor(const string &datafile, const string &outputfile)
{
    FILE *gp = popen("gnuplot -persistent", "w");

    fprintf(gp, "set terminal pdfcairo size 6,6 enhanced font 'CMU Serif,14' background '#ffffff'\n");
    fprintf(gp, "set output '%s'\n", outputfile.c_str());

    fprintf(gp, "set xlabel 'x[0](t)' font 'CMU Serif,30' offset 0,-1.\n");
    fprintf(gp, "set ylabel 'x[N](t)' font 'CMU Serif,30' offset -1.5,0\n");

    fprintf(gp, "set xtics 0.5 font 'CMU Serif,25' nomirror out\n");
    fprintf(gp, "set ytics 0.5 font 'CMU Serif,25' nomirror out\n");
    fprintf(gp, "set xtics format '%%.1f'\n");
    fprintf(gp, "set ytics format '%%.1f'\n");
    fprintf(gp, "set tics scale 1.0, 0.5\n");

    fprintf(gp, "set lmargin 12\nset rmargin 4\nset tmargin 2\nset bmargin 5\n");
    fprintf(gp, "set border 3 lw 1.5 lc rgb '#222222'\n");
    fprintf(gp, "set size ratio 1\n");

    fprintf(gp, "set datafile separator ','\n");

    fprintf(gp, "set grid xtics ytics lc rgb '#393939' lw 2 dt 1\n");
    fprintf(gp, "set grid back\n");

    fprintf(gp, "plot '%s' every ::1 using 1:2 with lines lw 0.001 lc rgb '#7a32a4' notitle\n", datafile.c_str());

    fflush(gp);
    pclose(gp);
}