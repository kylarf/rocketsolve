#ifndef ISENTROPIC_H
#define ISENTROPIC_H

#include <math.h>

typedef struct
{
    double A_Astar, gm;
} area_mach_par;

double area_mach(double M, area_mach_par par);

double p0_from_p(double p, double M, double gamma);

double T0_from_T(double T, double M, double gamma);

double rho0_from_rho(double rho, double M, double gamma);

double p_from_p0(double p0, double M, double gamma);

double T_from_T0(double T0, double M, double gamma);

double rho_from_rho0(double rho0, double M, double gamma);

#endif // ISENTROPIC_H
