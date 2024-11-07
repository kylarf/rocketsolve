#ifndef NOZZLE_H
#define NOZZLE_H

#include <math.h>

double calc_C_F(double p0, double p_e, double p_a, double Ae_At, double gm);

double calc_cstar(double T0, double Rspec, double gm);

#endif // NOZZLE_H
