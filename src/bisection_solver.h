#include <math.h>

#ifndef BISECTION_SOLVER_H
#define BISECTION_SOLVER_H
static inline int sign(double num)
{
    return (0.0 < num) - (num < 0.0);
}
#endif // BISECTION_SOLVER_H

#define SOLVER_CONCAT_INNER(A, B) A ## B
#define SOLVER_CONCAT(A, B) SOLVER_CONCAT_INNER(A, B)

#ifndef SOLVER_FUNC
#error SOLVER_FUNC must be defined before including "bisection_solver.h"
#endif

#ifndef SOLVER_NAME
    #define SOLVER_NAME SOLVER_CONCAT(bisect_solve_, SOLVER_FUNC)
#endif

#ifdef SOLVER_SIGNATURE
    #define COMMA_SIGNATURE , SOLVER_SIGNATURE
#else
    #define COMMA_SIGNATURE
#endif

#ifdef SOLVER_ARGS
    #define COMMA_ARGS , SOLVER_ARGS
#else
    #define COMMA_ARGS
#endif


double SOLVER_NAME(double a, double b, double tol, int nmax COMMA_SIGNATURE)
{
    int n = 1;
    double f_a, f_c, c;
    f_a = SOLVER_FUNC(a COMMA_ARGS);
    if (sign(f_a) != sign(SOLVER_FUNC(b COMMA_ARGS))) {
        while (n <= nmax) {
            c = (a + b) / 2.0;
            if ((f_c = SOLVER_FUNC(c COMMA_ARGS)) == 0 || (b - a)/2.0 < tol) {
                return c;
            }
            ++n;
            if (sign(f_c) == sign(f_a)) {
                a = c;
                f_a = f_c;
            }
            else {
                b = c;
            }
        }
    }
    return NAN;
}

#undef SOLVER_FUNC
#undef SOLVER_SIGNATURE
#undef COMMA_SIGNATURE
#undef SOLVER_ARGS
#undef COMMA_ARGS
#undef SOLVER_NAME
#undef SOLVER_CONCAT
#undef SOLVER_CONCAT_INNER
