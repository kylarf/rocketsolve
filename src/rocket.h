#ifndef ROCKET_H
#define ROCKET_H

#include <stdio.h>
#include <stdlib.h>

#define NUM_POINTS 1501

typedef struct PerfParams
{
    double p0, T0, h0, p_e, T_e, h_e, Ae_At, L_div, M_e, mdot, cstar, C_F, Isp, r_t;
} PerfParams;

typedef struct FlowProperties
{
    double *A, *M, *u, *p, *T, *rho, *h;
} FlowProperties;

typedef struct RocketInputs
{
    size_t n_points;
    double *x_loc, *nozzle_rad, *A_ratio;
    double F_t_nom, p_a, p_c, T_c, gamma, Mbar, Rspec, a, b, c, R_wtd, R_wtu,
           R_i, theta_i, CR;
} RocketInputs;

typedef struct Rocket
{
    RocketInputs *inputs;
    FlowProperties *flow_props;
    PerfParams *perf_params;
} Rocket;

Rocket *Rocket_alloc();

int Rocket_init(Rocket *self, FILE *inputs, FILE *nozzle_geometry);

void Rocket_free(Rocket *self);

void Rocket_compute_flow(Rocket *self);

void Rocket_compute_perf(Rocket *self);

void Rocket_write_summary(Rocket *self, FILE *summary);

int Rocket_write_rawdata(Rocket *self, char *datadir);

void gnuplot_plot_array(FILE *gnuplot, char *datadir, char *qty_name, char *ylabel);

void gnuplot_plot_with_stagn(FILE *gnuplot, char *datadir, char *qty_name, char *ylabel,
                             double stag, double xmin, double xmax);

int Rocket_make_plots(Rocket *self, char *datadir);

#endif // ROCKET_H
