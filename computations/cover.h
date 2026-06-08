#include <capd/capdlib.h>
#include <vector>

using namespace std;
using namespace capd;

/// @brief Checks whether a set s = x0 + B * r satisfies s =pm=> s,
///        then plots the visualisation of the covering relation
/// @param iPm      Interval Poincaré map defined for the set s.
/// @param A        Coordinate basis for values of the iPm.
/// @param B        Coordinate basis for the set s.
/// @param r        The vector of intervals centred around zero for defining the sides of the set s.
/// @param x0       IVector that is the centre of the set s.
bool checkSelfCover(IPoincareMap IPm, const IMatrix& A, const IMatrix& B, const IVector& r, const IVector& x0);

/// @brief Plots the rectangles set and P_set for visualisation of the covering relation
/// @param set          The set that will be covered.
/// @param P_set        The set Poincare map value of the original set that covers that set.
/// @param set_left     The left side of interval in the unstable direction of the original set.
/// @param P_set_left   The Pooincare map value of the set_left.
/// @param set_right    The right side of interval in the unstable direction of the original set.
/// @param P_set_right  The Pooincare map value of the set_right.
void plotRectangles(const IVector &set, const IVector &P_set,
                    const IVector &set_left, const IVector &P_set_left,
                    const IVector &set_right, const IVector &P_set_right);
