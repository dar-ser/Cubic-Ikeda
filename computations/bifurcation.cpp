#include "bifurcation.h"

void plotBifurcationDiagram(DMap &CubicIkeda, DPoincareMap &pm, const DVector &x,
                            double aStart, double aEnd, double aFrequency,
                            int noSteps, double chaoticA)
{
    FILE *gp = popen("gnuplot -persistent", "w");

    fprintf(gp, "set terminal pngcairo size 1800,1200 font 'Helvetica,11'\n");
    fprintf(gp, "set output 'images/bifurcation.png'\n");

    fprintf(gp, "set xlabel 'a' offset 0,-0.5\n");
    fprintf(gp, "set ylabel 'x(t-{/Symbol t})' offset 1,0\n");

    fprintf(gp, "set lmargin 10\nset rmargin 3\nset tmargin 2\nset bmargin 4\n");

    fprintf(gp, "set xrange [%f:%f]\n", aStart, aEnd);
    fprintf(gp, "set format x '%%.2f'\n");
    fprintf(gp, "set format y '%%.1f'\n");
    fprintf(gp, "set tics out nomirror\n");
    fprintf(gp, "set border 3\n");

    // vertical line at chaotic parameter
    fprintf(gp,
            "set arrow from %f, graph 0 to %f, graph 1 "
            "nohead lc rgb '#cc0000' lw 1 dt 2\n", // dashed red
            chaoticA, chaoticA);
    fprintf(gp,
            "set label 'a = %.3f' at %f, graph 0.97 "
            "left font 'Helvetica,9' tc rgb '#cc0000' offset 0.4,0\n",
            chaoticA, chaoticA);

    // plots the rest of the points

    fprintf(gp, "set pointsize 0.01\n");
    fprintf(gp, "plot '-' with points pt 7 lc rgb '#A7000000' notitle\n");

    double aIncrease = (aEnd - aStart) / aFrequency;
    int N = x.dimension() - 1;

    for (double a = aStart; a <= aEnd; a += aIncrease)
    {
        DVector returnPoint(x);
        CubicIkeda.setParameters({a});
        cout << "plotting bifurcation diagram for a = " << a << endl;

        try
        {
            for (int k = 0; k < 3. / 4. * noSteps; k++)
                returnPoint = pm(returnPoint); // skipping plotting for the first iterations
            for (int k = 3. / 4. * noSteps; k < noSteps; k++)
            {
                returnPoint = pm(returnPoint);
                fprintf(gp, "%.6f %.6f\n", a, returnPoint[N]); // plotting last iterations
            }
        }
        catch (exception &e)
        {
            cout << "Exception at a=" << a << ": " << e.what() << "\n";
            returnPoint = x;
        }
    }

    fprintf(gp, "e\n");
    fflush(gp);
    pclose(gp);
}

void computeBifurcationDiagram(DMap &CubicIkeda, DPoincareMap &pm, const DVector &x,
                               double aStart, double aEnd, double aFrequency,
                               int noSteps, const string &filename)
{
    FILE *data = fopen(filename.c_str(), "w");

    double aIncrease = (aEnd - aStart) / aFrequency;
    int N = x.dimension() - 1;

    for (double a = aStart; a <= aEnd; a += aIncrease)
    {
        DVector returnPoint(x);
        CubicIkeda.setParameters({a});

        cout << "computing for a = " << a << endl;

        try
        {
            for (int k = 0; k < 3.0 / 4.0 * noSteps; k++)
                returnPoint = pm(returnPoint); // we don't save first iterations

            for (int k = 3.0 / 4.0 * noSteps; k < noSteps; k++)
            {
                returnPoint = pm(returnPoint);
                fprintf(data, "%.10f %.10f\n", a, returnPoint[N]); // saving the rest of iterations
            }
        }
        catch (exception &e)
        {
            cout << "Exception at a = " << a << ": " << e.what() << endl;
        }
    }

    fclose(data);
}


void plotBifurcationDiagram(const string &datafile, const string &outputfile,
                            double aStart, double aEnd, double chaoticA)
{
    FILE *gp = popen("gnuplot -persistent", "w");

    fprintf(gp, "set terminal pngcairo size 2400,1400 enhanced font 'CMU Serif,28' background '#ffffff'\n");
    fprintf(gp, "set output '%s'\n", outputfile.c_str());

    fprintf(gp, "set xlabel '{/CMU-Serif-Italic a}' font 'CMU Serif,44' offset 0,-0.5\n");
    fprintf(gp, "set ylabel '{/CMU-Serif-Italic x}({/CMU-Serif-Italic t} - {/Symbol t})' font 'CMU Serif,44' offset 0.15,0\n");

    fprintf(gp, "set grid xtics ytics lc rgb '#cccccc' lw 0.5 dt 1\n");
    fprintf(gp, "set grid back\n");

    fprintf(gp, "set xtics font 'CMU Serif,26' nomirror out\n");
    fprintf(gp, "set ytics font 'CMU Serif,26' nomirror out\n");
    fprintf(gp, "set xtics format '%%.2f'\n");
    fprintf(gp, "set ytics format '%%.1f'\n");
    fprintf(gp, "set mxtics 5\n");
    fprintf(gp, "set mytics 5\n");
    fprintf(gp, "set tics scale 1.2, 0.6\n");

    fprintf(gp, "set lmargin 14\nset rmargin 5\nset tmargin 3\nset bmargin 6\n");
    fprintf(gp, "set xrange [%f:%f]\n", aStart, aEnd);
    fprintf(gp, "set border 3 lw 2.0 lc rgb '#222222'\n");

    fprintf(gp,
        "set arrow 1 from %f, graph 0 to %f, graph 1 "
        "nohead lc rgb '#bb2222' lw 4.0 dt (8,4) front\n",
        chaoticA, chaoticA);

    fprintf(gp,
        "set label 1 '{/CMU-Serif-Italic a}_{c} = %.3f' "
        "at %f, graph 0.95 left font 'CMU Serif,34' "
        "tc rgb '#bb2222' offset 0.7,0\n",
        chaoticA, chaoticA);

    fprintf(gp, "plot '%s' with points pt 7 ps 0.2 lc rgb '#000000' notitle, \n", datafile.c_str());

    fflush(gp);
    pclose(gp);
}