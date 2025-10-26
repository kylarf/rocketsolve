#include "nozzle.h"
#include "isentropic.h"

#include <math.h>

double M2_norm_shock(double M1, double gamma)
{ 
    return sqrt((1 + (gamma - 1)/2 * M1*M1)
                / (gamma*M1*M1 - (gamma - 1)/2));
}

double Me_norm_shock_in(double gm, double p0, double p_b, double Ae_At)
{
    return sqrt(-1/(gm-1)
                + sqrt(1/pow(gm-1, 2)
                       + 2/(gm-1) * pow(2/(gm+1), (gm+1)/(gm-1))
                       * pow(p0/p_b/Ae_At, 2)));
}

double rho_norm_shock(double M1, double gamma)
{
    return (gamma + 1) * M1*M1 / (2 + (gamma - 1)*M1*M1);
}

double T_norm_shock(double M1, double gamma)
{
    return (1 + 2*gamma/(gamma+1) * (M1*M1 - 1))
           * ((2 + (gamma - 1)*M1*M1) / ((gamma+1)*M1*M1));
}

double p_norm_shock(double M1, double gamma)
{
    return (1 + 2*gamma/(gamma+1) * (M1*M1 - 1));
}

double p0_norm_shock(double M1, double gamma)
{
    double numer = pow(((gamma+1)/2 * M1*M1) / (1 + (gamma-1)/2 * M1*M1), gamma/(gamma-1));
    double denom = pow(2*gamma/(gamma+1) * M1*M1 - (gamma-1)/(gamma+1), 1/(gamma-1));
    return numer / denom;
}

double calc_C_F(double p0, double p_e, double p_a, double Ae_At, double gm)
{
    return sqrt(2*gm*gm / (gm-1) * pow(2/(gm+1), (gm+1)/(gm-1))
                * (1 - pow(p_e/p0, (gm-1)/gm)))
                + (p_e - p_a) / p0 * Ae_At;
}

double calc_C_F_corr(double p0, double p_e, double p_a, double Ae_At, double gm,
                     double lm_div)
{
    return lm_div * sqrt(2*gm*gm / (gm-1) * pow(2/(gm+1), (gm+1)/(gm-1))
                * (1 - pow(p_e/p0, (gm-1)/gm)))
                + (p_e - p_a) / p0 * Ae_At;
}

double calc_cstar(double T0, double Rspec, double gm)
{
    return sqrt(1/gm * pow((gm+1)/2, (gm+1)/(gm-1)) * Rspec * T0);
}
