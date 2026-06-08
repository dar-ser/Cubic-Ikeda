#include "cover.h"
#include "newton.h"

bool checkSelfCover(IPoincareMap IPm, const IMatrix& A, const IMatrix& B, const IVector& r, const IVector& x0){
    int N = x0.dimension()-1;
    
    // we need to check :
    // 1. whether points laying on the left and right bounds
    //    in the unstable direrction (assumed to be the 1. coordinate)
    //      go outside the set
    // 2. whether points laying in the stable direction
    //      go inside the set

    // to answear the first question we need to define new centre points x0_left, x0_right
    // and new rectangles r_left_inter, r_right_inter for sets:
    //  s_left = x0_left + B * r_left_inter
    //  s_right = x0_right + B * r_right_inter

    IVector r_left(N + 1), r_left_point(N + 1), r_left_inter(N + 1);
    r_left[0] = r_left_point[0] = r_left_inter[0] = 0;
    r_left[1] = r_left_point[1] = r[1].leftBound();
    r_left_inter[1] = 0;
    
    IVector r_right(N + 1), r_right_point(N + 1), r_right_inter(N + 1);
    r_right[0] = r_right_point[0] = r_right_inter[0] = 0;
    r_right[1] = r_right_point[1] = r[1].rightBound();
    r_right_inter[1] = 0;
    
    for (int i = 2; i < N + 1; i++)
    {
        r_left_inter[i] = r_left[i] = r[i];
        r_right_inter[i] = r_right[i] = r[i];
        r_left_point[i] = 0;
        r_right_point[i] = 0;
    }

    IVector x0_left = x0 + B * r_left_point;
    
    x0_left[0] = 0;
    interval IReturnTime = 0;
    
    C0Rect2Set s_left(x0_left, B, r_left_inter);
    IVector P0_left = IPm(s_left, x0, A, IReturnTime);
    
    
    IVector x0_right = x0 + B * r_right_point;
    
    x0_right[0] = 0;
    IReturnTime = 0;
    
    C0Rect2Set s_right(x0_right, B, r_right_inter);
    IVector P0_right = IPm(s_right, x0, A, IReturnTime);   


    // checking the covering:

    bool leftCovering = P0_left[1].rightBound() < r[1].leftBound() ? true : false;
    bool rightCovering = r[1].rightBound() < P0_right[1].leftBound() ? true : false;
    cout << "Left covering: " << P0_left[1].rightBound() << " < " << r[1].leftBound() << " : " << leftCovering << endl;
    cout << "Right covering: " << r[1].rightBound() << " < " << P0_right[1].leftBound() << " : " << rightCovering << endl;
    
    C0Rect2Set s(x0, B, r);
    IVector P0 = IPm(s, x0, A, IReturnTime);
    
    
    
    // to answear the second question this is enough:
    
    bool stableCovering = true;
    
    for(int i = 2; i < r.dimension(); i++){
        if (!( r[i].leftBound() < P0[i].leftBound() && P0[i].rightBound() < r[i].rightBound() )){
            stableCovering = false;
            break;
        }
    }
    
    plotRectangles(r, P0, r_left, P0_left, r_right, P0_right);
    return leftCovering && rightCovering && stableCovering;
}



void plotRectangles(
    const IVector &set, const IVector &P_set,
    const IVector &set_left, const IVector &P_set_left,
    const IVector &set_right, const IVector &P_set_right)
{
    int N = set.dimension() - 1;

    FILE *gp = popen("gnuplot", "w");
    if (!gp)
        return;

    const int cols = 2;
    const int rows = (N - 1 + cols - 1) / cols;

    fprintf(gp,
        "set terminal pdfcairo size 16,10 enhanced font 'Helvetica,20' linewidth 2\n"
        "set output 'images/cover.pdf'\n"
        "set multiplot layout %d,%d rowsfirst spacing 0.06,0.06 margins 0.10,0.98,0.08,0.95\n",
        rows, cols);

    for (int k = 2; k <= N; ++k)
    {
        double xmin =
            std::min(
                std::min(set[1].leftBound(), P_set[1].leftBound()),
                std::min(
                    std::min(set_left[1].leftBound(), P_set_left[1].leftBound()),
                    std::min(set_right[1].leftBound(), P_set_right[1].leftBound())
                )
            );

        double xmax =
            std::max(
                std::max(set[1].rightBound(), P_set[1].rightBound()),
                std::max(
                    std::max(set_left[1].rightBound(), P_set_left[1].rightBound()),
                    std::max(set_right[1].rightBound(), P_set_right[1].rightBound())
                )
            );

        double ymin =
            std::min(
                std::min(set[k].leftBound(), P_set[k].leftBound()),
                std::min(
                    std::min(set_left[k].leftBound(), P_set_left[k].leftBound()),
                    std::min(set_right[k].leftBound(), P_set_right[k].leftBound())
                )
            );

        double ymax =
            std::max(
                std::max(set[k].rightBound(), P_set[k].rightBound()),
                std::max(
                    std::max(set_left[k].rightBound(), P_set_left[k].rightBound()),
                    std::max(set_right[k].rightBound(), P_set_right[k].rightBound())
                )
            );

        double xpad = 0.03 * (xmax - xmin);
        double ypad = 0.03 * (ymax - ymin);

        if (xpad == 0.0) xpad = 1e-12;
        if (ypad == 0.0) ypad = 1e-12;

        fprintf(gp,
            "set border lw 1.2\n"
            "set grid xtics ytics lw 0.5 dt 2\n"
            "set tics out nomirror scale 0.8\n"
            "set ylabel 'coord %d' font ',28' offset 0,0\n",
            k);

        fprintf(gp,
            "set xrange [%e:%e]\n"
            "set yrange [%e:%e]\n",
            xmin - xpad, xmax + xpad,
            ymin - ypad, ymax + ypad);

        // Remove tick labels (but keep tick marks/grid structure)
        fprintf(gp,
            "set format x ''\n"
            "set format y ''\n");

        fprintf(gp,
            "plot "
            "'-' with filledcurves closed lc rgb '#D55E00' fs solid 0.30 notitle, "
            "'-' with filledcurves closed lc rgb '#67dd5a' fs solid 0.30 notitle, "
            "'-' with lines lc rgb '#5676e9' lw 2.5 notitle, "
            "'-' with lines lc rgb '#0900b2' lw 2.5 notitle, "
            "'-' with lines lc rgb '#ff42a1' lw 2.5 notitle, "
            "'-' with lines lc rgb '#882255' lw 2.5 notitle\n");

        fprintf(gp,
            "%e %e\n%e %e\n%e %e\n%e %e\n%e %e\ne\n",
            set[1].leftBound(),  set[k].leftBound(),
            set[1].rightBound(), set[k].leftBound(),
            set[1].rightBound(), set[k].rightBound(),
            set[1].leftBound(),  set[k].rightBound(),
            set[1].leftBound(),  set[k].leftBound());

        fprintf(gp,
            "%e %e\n%e %e\n%e %e\n%e %e\n%e %e\ne\n",
            P_set[1].leftBound(),  P_set[k].leftBound(),
            P_set[1].rightBound(), P_set[k].leftBound(),
            P_set[1].rightBound(), P_set[k].rightBound(),
            P_set[1].leftBound(),  P_set[k].rightBound(),
            P_set[1].leftBound(),  P_set[k].leftBound());

        fprintf(gp,
            "%e %e\n%e %e\ne\n",
            set_left[1].leftBound(),  set_left[k].leftBound(),
            set_left[1].rightBound(), set_left[k].rightBound());

        fprintf(gp,
            "%e %e\n%e %e\ne\n",
            P_set_left[1].leftBound(),  P_set_left[k].leftBound(),
            P_set_left[1].rightBound(), P_set_left[k].rightBound());

        fprintf(gp,
            "%e %e\n%e %e\ne\n",
            set_right[1].leftBound(),  set_right[k].leftBound(),
            set_right[1].rightBound(), set_right[k].rightBound());

        fprintf(gp,
            "%e %e\n%e %e\ne\n",
            P_set_right[1].leftBound(),  P_set_right[k].leftBound(),
            P_set_right[1].rightBound(), P_set_right[k].rightBound());
    }

    fprintf(gp, "unset multiplot\n");
    fprintf(gp, "unset output\n");

    fflush(gp);
    pclose(gp);
}