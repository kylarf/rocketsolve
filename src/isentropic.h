#ifndef ISENTROPIC_H
#define ISENTROPIC_H

#include <math.h>

double area_mach(double M, double A_Astar, double gm);

double a_astar(double M, double gm);

double p0_from_p(double p, double M, double gamma);

double T0_from_T(double T, double M, double gamma);

double rho0_from_rho(double rho, double M, double gamma);

double p_from_p0(double p0, double M, double gamma);

double T_from_T0(double T0, double M, double gamma);

double rho_from_rho0(double rho0, double M, double gamma);

double pran_mey(double M, double gamma);

#endif // ISENTROPIC_H
