#include <math.h>

double calc_C_F(double p0, double p_e, double p_a, double Ae_At, double gm)
{
    return sqrt(2*gm*gm / (gm-1) * pow(2/(gm+1), (gm+1)/(gm-1))
                * (1 - pow(p_e/p0, (gm-1)/gm)))
                + (p_e - p_a) / p0 * Ae_At;
}

double calc_cstar(double T0, double Rspec, double gm)
{
    return sqrt(1/gm * pow((gm+1)/2, (gm+1)/(gm-1)) * Rspec * T0);
}
