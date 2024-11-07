#include <malloc/_malloc.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "rocket.h"
#include "fileio.h"
#include "solver.h"
#include "isentropic.h"
#include "nozzle.h"
#include "constants.h"

int read_inputs(FILE *input_file, RocketInputs *rocket_inputs)
{
    // Read a parameter from input file
    #define READ_INPUT(field) \
    if (fgets(line, sizeof(line), input_file) != NULL) { \
        temp = strtod(line, &endptr); \
        if (line == endptr) { \
            return 1; \
        } \
        else { \
            rocket_inputs->field = temp; \
        } \
    } \
    else { \
        return 1; \
    }
    double temp;
    char line[256];
    char *endptr;

    // ensure we are at beginning of file
    fseek(input_file, 0, SEEK_SET);

    // Read parameters. Assumes fixed order of input parameters in file
    READ_INPUT(F_t_nom);
    READ_INPUT(p_a);
    READ_INPUT(p_c);
    READ_INPUT(T_c);
    READ_INPUT(gamma);
    READ_INPUT(Mbar);
    READ_INPUT(a);
    READ_INPUT(b);
    READ_INPUT(c);
    READ_INPUT(R_wtd);
    READ_INPUT(R_wtu);
    READ_INPUT(R_i);
    READ_INPUT(theta_i);
    READ_INPUT(CR);

    return 0;
    #undef READ_INPUT
}

Rocket *Rocket_alloc()
{
    Rocket *rocket = malloc(sizeof(*rocket));
    rocket->inputs = calloc(1, sizeof(RocketInputs));
    rocket->flow_props = calloc(1, sizeof(FlowProperties));
    rocket->perf_params = calloc(1, sizeof(PerfParams));
    return rocket;
}

int Rocket_init(Rocket *self, FILE *inputs, FILE *nozzle_geometry)
{
    RocketInputs *rocket_inputs = self->inputs;
    FlowProperties *flow_props = self->flow_props;
    
    // count lines in geometry file == number of eval points + 1
    size_t n_points = count_lines(nozzle_geometry) - 1;
    rocket_inputs->n_points = n_points;

    // read values from input file
    if (read_inputs(inputs, rocket_inputs) != 0) {
        return 1;
    }
    rocket_inputs->Rspec = R_BAR / rocket_inputs->Mbar;

    // allocate memory for flow property arrays
    rocket_inputs->x_loc = calloc(n_points, sizeof(double));
    rocket_inputs->nozzle_rad = calloc(n_points, sizeof(double));
    rocket_inputs->A_ratio = calloc(n_points, sizeof(double));
    flow_props->A = calloc(n_points, sizeof(double));
    flow_props->M = calloc(n_points, sizeof(double));
    flow_props->u = calloc(n_points, sizeof(double));
    flow_props->p = calloc(n_points, sizeof(double));
    flow_props->T = calloc(n_points, sizeof(double));
    flow_props->rho = calloc(n_points, sizeof(double));
    flow_props->h = calloc(n_points, sizeof(double));

    // read nozzle geometry from file into x and r arrays
    read_geometry(nozzle_geometry, n_points, rocket_inputs->x_loc, rocket_inputs->nozzle_rad);

    // compute and populate area ratio array from radii
    for (size_t i = 0; i < n_points; ++i) {
        double r_i = rocket_inputs->nozzle_rad[i];
        rocket_inputs->A_ratio[i] = r_i * r_i;
    }
    return 0;
}

void Rocket_free(Rocket *self)
{
    RocketInputs *inputs = self->inputs;
    FlowProperties *flow_props = self->flow_props;

    free(inputs->x_loc);
    free(inputs->nozzle_rad);
    free(inputs->A_ratio);
    free(inputs);

    free(flow_props->A);
    free(flow_props->M);
    free(flow_props->u);
    free(flow_props->p);
    free(flow_props->T);
    free(flow_props->rho);
    free(flow_props->h);
    free(flow_props);

    free(self->perf_params);
    free(self);
}

IMPL_SOLVER(area_mach, area_mach_par, par);

void Rocket_compute_flow(Rocket *self)
{
    // Assumes that Rocket_init() has been called. UNDEFINED BEHAVIOR IF NOT as
    // arrays are unallocated. Will probably (hopefully) segfault.
    RocketInputs *inputs = self->inputs;
    FlowProperties *flow_props = self->flow_props;
    PerfParams *perf_params = self->perf_params;

    size_t n_points = inputs->n_points;
    double rho0, c_local, a, b, gamma = inputs->gamma;
    double c_p = gamma * inputs->Rspec / (gamma - 1);

    // initial pass to compute Mach numbers
    for (size_t i = 0; i < n_points; ++i) {
        // choose subsonic or supersonic based on position relative to throat
        if (inputs->x_loc[i] < 0.0) {
            a = 0;
            b = 1;
        }
        else {
            a = 1;
            b = 20;
        }
        flow_props->M[i] = solve_area_mach(a, b, 1e-8, 50,
                                           (area_mach_par){inputs->A_ratio[i], gamma});
    }
    perf_params->p0 = p0_from_p(inputs->p_c, flow_props->M[0], gamma);
    perf_params->T0 = T0_from_T(inputs->T_c, flow_props->M[0], gamma);
    perf_params->h0 = c_p * perf_params->T0;
    rho0 = perf_params->p0 / inputs->Rspec / perf_params->T0;

    // finally calculate flow properties through nozzle
    for (size_t i = 0; i < n_points; ++i) {
        // flow property calculations
        flow_props->p[i] = p_from_p0(perf_params->p0, flow_props->M[i], gamma);
        flow_props->T[i] = T_from_T0(perf_params->T0, flow_props->M[i], gamma);
        flow_props->rho[i] = rho_from_rho0(rho0, flow_props->M[i], gamma);
        flow_props->h[i] = c_p * flow_props->T[i];
        // local speed of sound
        c_local = sqrt(gamma * inputs->Rspec * flow_props->T[i]);
        flow_props->u[i] = c_local * flow_props->M[i];
    }
}

void Rocket_compute_perf(Rocket *self)
{
    RocketInputs *inputs = self->inputs;
    FlowProperties *flow_props = self->flow_props;
    PerfParams *perf_params = self->perf_params;
    size_t n_points = inputs->n_points;
    double Astar;

    // exit properties
    perf_params->p_e = flow_props->p[n_points-1];
    perf_params->T_e = flow_props->T[n_points-1];
    perf_params->h_e = flow_props->h[n_points-1];
    perf_params->Ae_At = inputs->A_ratio[n_points-1];
    perf_params->M_e = flow_props->M[n_points-1];

    // performance parameters
    perf_params->C_F = calc_C_F(perf_params->p0, perf_params->p_e, inputs->p_a,
                                perf_params->Ae_At, inputs->gamma);
    perf_params->cstar = calc_cstar(perf_params->T0, inputs->Rspec, inputs->gamma);
    perf_params->Isp = perf_params->cstar * perf_params->C_F / G_0;
    
    // absolute nozzle throat area and radius, nozzle diverging length
    Astar = inputs->F_t_nom / perf_params->p0 / perf_params->C_F;
    perf_params->r_t = sqrt(Astar * M_1_PI);
    perf_params->L_div = inputs->x_loc[n_points-1] * perf_params->r_t;
    for (size_t i = 0; i < n_points; ++i) {
        flow_props->A[i] = Astar * inputs->A_ratio[i];
    }

    perf_params->mdot = perf_params->p0 * Astar / perf_params->cstar;
}

void Rocket_write_summary(Rocket *self, FILE *summary)
{
    PerfParams *perf_params = self->perf_params;
    fprintf(summary, "p0 = %#.6g Pa\n", perf_params->p0);
    fprintf(summary, "T0 = %#.6g K\n", perf_params->T0);
    fprintf(summary, "h0 = %#.6g J/kg\n", perf_params->h0);
    fprintf(summary, "pe = %#.6g Pa\n", perf_params->p_e);
    fprintf(summary, "Te = %#.6g K\n", perf_params->T_e);
    fprintf(summary, "he = %#.6g J/kg\n", perf_params->h_e);
    fprintf(summary, "Ae/At = %#.6g\n", perf_params->Ae_At);
    fprintf(summary, "Ld = %#.6g m\n", perf_params->L_div);
    fprintf(summary, "Me = %#.6g\n", perf_params->M_e);
    fprintf(summary, "mdot = %#.6g kg/s\n", perf_params->mdot);
    fprintf(summary, "cstar = %#.6g m/s\n", perf_params->cstar);
    fprintf(summary, "C_F = %#.6g\n", perf_params->C_F);
    fprintf(summary, "Isp = %#.6g s\n", perf_params->Isp);
}

int cat_path_write(char *path_buf, char *datadir, char *filename, double *arr, int N)
{
    // make buffer appear to be empty (i.e. ending with null terminator)
    path_buf[0] = '\0';
    // construct path to data file
    strcat(path_buf, datadir);
    strcat(path_buf, "/");
    strcat(path_buf, filename);
    return write_array(path_buf, arr, N);
}

int Rocket_write_rawdata(Rocket *self, char *datadir)
{
    int status = 0;
    char path_buf[4096];
    RocketInputs *inputs = self->inputs;
    FlowProperties *flow_props = self->flow_props;
    size_t n_points = self->inputs->n_points;

    status |= cat_path_write(path_buf, datadir, "area.dat", flow_props->A, n_points);
    status |= cat_path_write(path_buf, datadir, "flow_velocity.dat", flow_props->u,
                             n_points);
    status |= cat_path_write(path_buf, datadir, "mach.dat", flow_props->M, n_points);
    status |= cat_path_write(path_buf, datadir, "specific_enthalpy.dat",
                             flow_props->h, n_points);
    status |= cat_path_write(path_buf, datadir, "static_density.dat",
                             flow_props->rho, n_points);
    status |= cat_path_write(path_buf, datadir, "static_pressure.dat",
                             flow_props->p, n_points);
    status |= cat_path_write(path_buf, datadir, "static_temperature.dat",
                             flow_props->T, n_points);
    status |= cat_path_write(path_buf, datadir, "x_position.dat",
                             inputs->x_loc, n_points);
    return status;
}

void gnuplot_plot_array(FILE *gnuplot, char *datadir, char *qty_name, char *ylabel)
{
    fputs("set terminal svg\n", gnuplot);
    fprintf(gnuplot, "set output '%s/%s.svg'\n", datadir, qty_name);
    fputs("set size ratio 0.5\n", gnuplot);
    fprintf(gnuplot, "set ylabel '%s'\n", ylabel);
    fprintf(gnuplot,
            "plot \"< paste %s/x_position.dat %s/%s.dat\" using 1:2 with lines notitle\n",
            datadir, datadir, qty_name);
}

void gnuplot_plot_with_stagn(FILE *gnuplot, char *datadir, char *qty_name, char *ylabel,
                             double stag, double xmin, double xmax)
{
    fputs("set terminal svg\n", gnuplot);
    fprintf(gnuplot, "set output '%s/%s.svg'\n", datadir, qty_name);
    fputs("set size ratio 0.5\n", gnuplot);
    fprintf(gnuplot, "set ylabel '%s'\n", ylabel);
    fprintf(gnuplot,
            "plot \"< paste %s/x_position.dat %s/%s.dat\" using 1:2 with lines notitle, "
            "'-' using 1:2 with lines title 'stagnation property'\n",
            datadir, datadir, qty_name);
    fprintf(gnuplot, "%.10g %.10g\n", xmin, stag);
    fprintf(gnuplot, "%.10g %.10g\ne\n", xmax, stag);
}

int Rocket_make_plots(Rocket *self, char *datadir)
{
    int status = 1;
    FILE *gnuplot = popen("gnuplot", "w");
    if (gnuplot == NULL) {
        goto gnuplot_failed;
    }
    
    gnuplot_plot_array(gnuplot, datadir, "mach", "Mach number M");
    gnuplot_plot_array(gnuplot, datadir, "area", "Area A (m^2)");
    gnuplot_plot_array(gnuplot, datadir, "flow_velocity", "Flow velocity u (m/s)");
    gnuplot_plot_array(gnuplot, datadir, "static_density",
                       "Static density rho (kg/m^3)");
    gnuplot_plot_with_stagn(gnuplot, datadir, "static_pressure", "Static pressure p (Pa)",
                            self->perf_params->p0, -5.0, 10.0);
    gnuplot_plot_with_stagn(gnuplot, datadir, "static_temperature",
                            "Static temperature T (K)",
                            self->perf_params->T0, -5.0, 10.0);
    gnuplot_plot_with_stagn(gnuplot, datadir, "specific_enthalpy",
                            "Specific enthalpy h (J/kg)",
                            self->perf_params->h0, -5.0, 10.0);

    pclose(gnuplot);
    status = 0;

gnuplot_failed:
    return status;
}
