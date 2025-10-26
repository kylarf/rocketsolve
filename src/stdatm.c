#include "stdatm.h"

#include <math.h>

static const double R_EARTH = 6369.0;   // radius of Earth (km)
static const double GMR = 34.163195;    // hydrostatic constant
static const int N_TAB = 8;             // length of tables
static const double P_SEA = 101325;     // pressure at sea level (Pa)
static const double T_SEA = 288.15;     // temperature at sea level (K)
static const double RHO_SEA = 1.2250;   // density at sea level (kg m^-3)

static const double h_tab[] = {0.0, 11.0, 20.0, 32.0, 47.0, 51.0, 71.0, 84.852};
static const double t_tab[] = {288.15, 216.65, 216.65, 228.65, 270.65, 270.65,
                               214.65, 186.946};
static const double p_tab[] = {1.0, 2.2336110E-1, 5.4032950E-2, 8.5666784E-3,
                               1.0945601E-3, 6.6063531E-4, 3.9046834E-5, 3.68501E-6};
static const double g_tab[] = {-6.5, 0.0, 1.0, 2.8, 0, -2.8, -2.0, 0.0};


void stdatm(double alt, double *p, double *T, double *rho)
{
    double h = alt * R_EARTH / (alt + R_EARTH); // geometric to geopotential altitude

    int i = 0, j = N_TAB;
    while (j > i + 1) {
        int k = (i + j) / 2;
        if (h < h_tab[k]) {
            j = k;
        }
        else {
            i = k;
        }
    }

    double t_grad = g_tab[i];                   // temp gradiant of local layer
    double t_base = t_tab[i];                   // base temp of local layer
    double delta_h = h - h_tab[i];              // height above local base
    double t_local = t_base + t_grad * delta_h; // local temp
    double theta = t_local / t_tab[0];          // temp ratio
    
    double delta, sigma;
    if (0.0 == t_grad) {
        delta = p_tab[i] * exp(-GMR * delta_h / t_base);
    }
    else {
        delta = p_tab[i] * pow(t_base / t_local, GMR / t_grad);
    }
    sigma = delta / theta;

    // convert from normalized values to SI units and set return vals
    *p = delta * P_SEA;
    *T = theta * T_SEA;
    *rho = sigma * RHO_SEA;
}
