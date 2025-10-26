#include <math.h>

#include "isentropic.h"

double area_mach(double M, area_mach_par par)
{
    return 1/M * pow(( 2/(par.gm+1) * (1 + (par.gm-1)/2 * M*M) ), (par.gm+1)/2/(par.gm-1))
           - par.A_Astar;
}

double a_astar(double M, double gm)
{
    return 1/M * pow(( 2/(gm+1) * (1 + (gm-1)/2 * M*M) ), (gm+1)/2/(gm-1));
}

double p0_from_p(double p, double M, double gamma)
{
    return p * pow(1 + (gamma-1)/2 * M*M, gamma/(gamma-1));
}

double T0_from_T(double T, double M, double gamma)
{
    return T * (1 + (gamma-1)/2 * M*M);
}

double rho0_from_rho(double rho, double M, double gamma)
{
    return rho * pow(1 + (gamma-1)/2 * M*M, 1/(gamma-1));
}

double p_from_p0(double p0, double M, double gamma)
{
    return p0 / pow(1 + (gamma-1)/2 * M*M, gamma/(gamma-1));
}

double T_from_T0(double T0, double M, double gamma)
{
    return T0 / (1 + (gamma-1)/2 * M*M);
}

double rho_from_rho0(double rho0, double M, double gamma)
{
    return rho0 / pow(1 + (gamma-1)/2 * M*M, 1/(gamma-1));
}

double pran_mey(double M, double gm)
{
    return sqrt((gm+1)/(gm-1)) * atan(sqrt((gm-1)/(gm+1) * (M*M-1))) - atan(sqrt(M*M-1));
}
