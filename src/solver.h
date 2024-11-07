#ifndef SOLVER_H
#define SOLVER_H

#include <math.h>

inline int sgn(double num)
{
    return (0.0 < num) - (num < 0.0);
}

#define DECL_0() 
#define DECL_2(t1, v1) , t1 v1
#define DECL_N(_2,_1,_0,N,...) DECL##N
#define DECL(...) DECL_N(_0,##__VA_ARGS__,_2,_1,_0)(__VA_ARGS__)

#define CALL_0()
#define CALL_2(t1, v1) , v1        
#define CALL_N(_2,_1,_0,N,...) CALL##N
#define CALL(...)  CALL_N(_0,##__VA_ARGS__,_2,_1,_0)(__VA_ARGS__)

#define IMPL_SOLVER(_f, ...) \
double solve_##_f(double _a, double _b, double _tol, int _nmax DECL(__VA_ARGS__)) \
{ \
    int _n = 1; \
    double _f_a, _f_c, _c; \
    _f_a = _f(_a CALL(__VA_ARGS__)); \
    if (sgn(_f_a) != sgn(_f(_b CALL(__VA_ARGS__)))) { \
        while (_n <= _nmax) { \
            _c = (_a + _b) / 2.0; \
            if ((_f_c = _f(_c CALL(__VA_ARGS__))) == 0 || (_b - _a)/2.0 < _tol) { \
                return _c; \
            } \
            ++_n; \
            if (sgn(_f_c) == sgn(_f_a)) { \
                _a = _c; \
                _f_a = _f_c; \
            } \
            else { \
                _b = _c; \
            } \
        } \
    } \
    return NAN; \
}

#endif // SOLVER_H
