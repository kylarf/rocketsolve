#include <stdio.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <string.h>
#include <limits.h>

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
    //rocket->inputs = calloc(1, sizeof(RocketInputs));
    //rocket->flow_props = calloc(1, sizeof(FlowProperties));
    //rocket->perf_params = calloc(1, sizeof(PerfParams));
    return rocket;
}

int Rocket_init(Rocket *self, FILE *inputs, FILE *nozzle_geometry)
{
    RocketInputs *rocket_inputs = &self->inputs;
    FlowProperties *flow_props = &self->flow_props;

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
    // also find throat index
    size_t i_th = 0;
    double r_i, r_min = rocket_inputs->nozzle_rad[0];
    for (size_t i = 0; i < n_points; ++i) {
        r_i = rocket_inputs->nozzle_rad[i];
        if (r_i < r_min) {
            r_min = r_i;
            i_th = i;
        }
        rocket_inputs->A_ratio[i] = r_i * r_i;
    }
    rocket_inputs->throat_idx = i_th;

    return 0;
}

Rocket *Rocket_from_inputs(RocketInputs *inputs)
{
    Rocket *rocket = Rocket_alloc();
    size_t n_points = inputs->n_points;
    size_t n_bytes = n_points * sizeof(double);

    // copy input parameter struct to self
    rocket->inputs = *inputs;

    // allocate memory and initialize nozzle geometry arrays
    rocket->inputs.x_loc = calloc(n_points, sizeof(double));
    rocket->inputs.nozzle_rad = calloc(n_points, sizeof(double));
    rocket->inputs.A_ratio = calloc(n_points, sizeof(double));
    memcpy(rocket->inputs.x_loc, inputs->x_loc, n_bytes);
    memcpy(rocket->inputs.nozzle_rad, inputs->nozzle_rad, n_bytes);
    memcpy(rocket->inputs.A_ratio, inputs->A_ratio, n_bytes);

    // allocate memory for flow property arrays
    rocket->flow_props.A = calloc(n_points, sizeof(double));
    rocket->flow_props.M = calloc(n_points, sizeof(double));
    rocket->flow_props.u = calloc(n_points, sizeof(double));
    rocket->flow_props.p = calloc(n_points, sizeof(double));
    rocket->flow_props.T = calloc(n_points, sizeof(double));
    rocket->flow_props.rho = calloc(n_points, sizeof(double));
    rocket->flow_props.h = calloc(n_points, sizeof(double));

    return rocket;
}

void Rocket_free(Rocket *self)
{
    RocketInputs *inputs = &self->inputs;
    FlowProperties *flow_props = &self->flow_props;

    free(inputs->x_loc);
    free(inputs->nozzle_rad);
    free(inputs->A_ratio);

    free(flow_props->A);
    free(flow_props->M);
    free(flow_props->u);
    free(flow_props->p);
    free(flow_props->T);
    free(flow_props->rho);
    free(flow_props->h);
    free(self);
}

IMPL_SOLVER(area_mach, area_mach_par, par);

void Rocket_compute_flow(Rocket *self)
{
    // Assumes that Rocket_init() has been called. UNDEFINED BEHAVIOR IF NOT as
    // arrays are unallocated. Will probably (hopefully) segfault.
    RocketInputs *inputs = &self->inputs;
    FlowProperties *flow_props = &self->flow_props;
    PerfParams *perf_params = &self->perf_params;

    size_t n_points = inputs->n_points;
    double rho0, c_local, a, b, gamma = inputs->gamma, p_b = inputs->p_a;
    double c_p = gamma * inputs->Rspec / (gamma - 1);
    double Me_sub, Me_sup, pe_sub, pe_sup;

    // initial pass to compute Mach numbers
    for (size_t i = 0; i < n_points; ++i) {
        // choose subsonic or supersonic based on whether flow has choked
        // (achieved sonic velocity)
        if (i >= inputs->throat_idx) {
            a = 1;
            b = 20; // something is probably wrong if Mach number is this high
        }
        else {
            a = 0;
            b = 1;
        }
        flow_props->M[i] = solve_area_mach(a, b, EPS, 50,
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

    handle_shocks(self);
}

void handle_shocks(Rocket *self)
{
    RocketInputs *inputs = &self->inputs;
    FlowProperties *flow_props = &self->flow_props;
    PerfParams *perf_params = &self->perf_params;
    
    size_t n_points = inputs->n_points;
    double Me_sub, Me_sup, pe_sub, pe_sup, pe2;
    double gamma = inputs->gamma, p_b = inputs->p_a;
    double c_p = gamma * inputs->Rspec / (gamma - 1);

    // back-pressure comparisons now
    Me_sub = solve_area_mach(0, 1, EPS, 50,
                             (area_mach_par){inputs->A_ratio[n_points-1], gamma});
    Me_sup = solve_area_mach(1, 20, EPS, 50, 
                             (area_mach_par){inputs->A_ratio[n_points-1], gamma});
    pe_sub = p_from_p0(perf_params->p0, Me_sub, gamma);
    pe_sup = p_from_p0(perf_params->p0, Me_sup, gamma);
    
    pe2 = pe_sup * p_norm_shock(Me_sup, gamma);

    printf("p_b,sub = %#.10g Pa\np_b,sup = %#.10g Pa\n", pe_sub, pe_sup);
    printf("p_b,exitshock = %#.10g Pa\n", pe2);
    
    if (fabs(p_b - pe_sub) < 1) {
        // isentropic choked subsonic expansion
        puts("Isentropic choked subsonic expansion.");
        flow_props->expansion.flavor = SONIC_THROAT;
        // initial pass to re-compute Mach numbers, subsonic
        for (size_t i = inputs->throat_idx; i < n_points; ++i) {
            flow_props->M[i] = solve_area_mach(0, 1, EPS, 50,
                                               (area_mach_par){inputs->A_ratio[i], gamma});
        }
        // recalculate flow properties through nozzle
        for (size_t i = inputs->throat_idx; i < n_points; ++i) {
            // flow property calculations
            flow_props->p[i] = p_from_p0(perf_params->p0, flow_props->M[i], gamma);
            flow_props->T[i] = T_from_T0(perf_params->T0, flow_props->M[i], gamma);
            double rho0 = perf_params->p0 / inputs->Rspec / perf_params->T0;
            flow_props->rho[i] = rho_from_rho0(rho0, flow_props->M[i], gamma);
            flow_props->h[i] = c_p * flow_props->T[i];
            // local speed of sound
            double c_local = sqrt(gamma * inputs->Rspec * flow_props->T[i]);
            flow_props->u[i] = c_local * flow_props->M[i];
        }
    }
    else if (fabs(p_b - pe_sup) < 1) {
        // isentropic supersonic expansion
        puts("Ideal isentropic expansion (supersonic).");
        flow_props->expansion.flavor = IDEAL;
    }
    else if (p_b > pe_sub) {
        // subsonic throughout
        puts("Subsonic nozzle flow throughout (unchoked).");
        flow_props->expansion.flavor = SUBSONIC;
    }
    else if (pe_sub > p_b && p_b > pe_sup) {
        double M2 = M2_norm_shock(Me_sup, gamma);

        if (fabs(pe2 - p_b) < 1) {
            // shock sitting at exit
            flow_props->expansion.flavor = NORM_SHK_EXIT;
            ExpNormShkExit *exp_props = &flow_props->expansion.props.norm_shk_exit;
            exp_props->shock_loc = inputs->x_loc[n_points-1];
            printf("Normal shock at exit, x_loc = %.2f\n", exp_props->shock_loc);
            flow_props->M[n_points-1] = M2;
            flow_props->p[n_points-1] = pe2;
            flow_props->rho[n_points-1] = flow_props->rho[n_points-1]
                                          * rho_norm_shock(Me_sup, gamma);
            flow_props->T[n_points-1] = flow_props->T[n_points-1]
                                        * T_norm_shock(Me_sup, gamma);
            flow_props->h[n_points-1] = c_p * flow_props->T[n_points-1];
            double c_exit = sqrt(gamma * inputs->Rspec * flow_props->T[n_points-1]);
            flow_props->u[n_points-1] = c_exit * flow_props->M[n_points-1];
        }
        else if (p_b > pe2) {
            // very overexpanded, normal shock inside
            flow_props->expansion.flavor = NORM_SHK_IN;
            ExpNormShkIn *exp_props = &flow_props->expansion.props.norm_shk_in;
            double M_e = Me_norm_shock_in(gamma, perf_params->p0, p_b,
                                          inputs->A_ratio[n_points-1]);
            //printf("M_e = %f\n", M_e);
            double Astar_ratio = a_astar(M_e, gamma) / inputs->A_ratio[n_points-1];
            //printf("Astar_ratio = %f\n", Astar_ratio);
            double p0_ratio, err_min = INFINITY, abs_err;
            size_t i_shk;
            for (size_t i = inputs->throat_idx; i < n_points; ++i) {
                p0_ratio = p0_norm_shock(flow_props->M[i], gamma);
                //printf("p0 ratio = %f\n", p0_ratio);
                abs_err = fabs(p0_ratio - Astar_ratio);
                if (abs_err < err_min) {
                    err_min = abs_err;
                    i_shk = i;
                }
            }
            exp_props->shock_loc = inputs->x_loc[i_shk];
            puts("Overexpanded; normal shock inside nozzle.");
            printf("Normal shock x_loc = %.2f\n", exp_props->shock_loc);
            double M1 = flow_props->M[i_shk];
            p0_ratio = p0_norm_shock(flow_props->M[i_shk], gamma);

            // recompute Mach numbers downstream of shock
            for (size_t i = i_shk; i < n_points; ++i) {
                // Mach number must be subsonic after shock
                flow_props->M[i] = solve_area_mach(
                    0, 1, EPS, 50, (area_mach_par){inputs->A_ratio[i]/Astar_ratio, gamma}
                );
            }

            double rho2 = flow_props->rho[i_shk] * rho_norm_shock(M1, gamma);
            double rho0 = rho0_from_rho(rho2, flow_props->M[i_shk], gamma);

            // finally recalculate flow properties through nozzle
            for (size_t i = i_shk; i < n_points; ++i) {
                // flow property calculations
                // use updated stagnation pressure
                flow_props->p[i] = p_from_p0(perf_params->p0*p0_ratio,
                                             flow_props->M[i], gamma);
                // stagnation temp does not change
                flow_props->T[i] = T_from_T0(perf_params->T0, flow_props->M[i],
                                             gamma);
                flow_props->rho[i] = rho_from_rho0(rho0, flow_props->M[i], gamma);
                flow_props->h[i] = c_p * flow_props->T[i];
                // local speed of sound
                double c_local = sqrt(gamma * inputs->Rspec * flow_props->T[i]);
                flow_props->u[i] = c_local * flow_props->M[i];
            }
        }
        else {
            // slightly overexpanded. no need to recompute properties, just
            // angle
            puts("Slightly overexpanded; oblique shock outside exit.");
            flow_props->expansion.flavor = OBL_SHK;
            double M_e = flow_props->M[n_points-1];
            double p_b_p_e = p_b / pe_sup;
            // solved analytically
            double beta_rad = asin(sqrt(1/M_e/M_e * ((gamma+1)/2/gamma * (p_b_p_e - 1) + 1)));
            // convert to degrees
            flow_props->expansion.props.obl_shk.shk_angle = beta_rad * 180 * M_1_PI;
            printf("Oblique shock angle = %#.6g deg\n",
                   flow_props->expansion.props.obl_shk.shk_angle);
        }
    }
    else if (pe_sup > p_b) {
        // underexpanded. no need to recompute inside nozzle though.
        puts("Underexpanded, expansion waves outside nozzle.");
        flow_props->expansion.flavor = EXP_FAN;
        double M_e = flow_props->M[n_points-1];
        double p0_pb = perf_params->p0 / p_b;
        double M_2 = sqrt(2/(gamma-1) * (pow(p0_pb, gamma/(gamma-1)) - 1));
        double turn_rad = pran_mey(M_2, gamma) - pran_mey(M_e, gamma);
        flow_props->expansion.props.expfan.turn_angle = turn_rad * 180 * M_1_PI;
        printf("Expansion wave turning angle = %#.6g deg\n",
               flow_props->expansion.props.expfan.turn_angle);
    }
    puts("");
}

void Rocket_compute_perf(Rocket *self)
{
    RocketInputs *inputs = &self->inputs;
    FlowProperties *flow_props = &self->flow_props;
    PerfParams *perf_params = &self->perf_params;
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

    double x_end = inputs->x_loc[n_points-1];
    double th_d = atan(inputs->a);
    double x_c0 = sin(th_d) * inputs->R_wtd;
    // derivative obtained analytically from functional form
    double dr_dx = inputs->a - 2*inputs->b*(x_end-x_c0) - 3*inputs->c*pow(x_end-x_c0,2);
    perf_params->th_e = atan(dr_dx) * 180 * M_1_PI;
    perf_params->lm_div = 0.5 * (1 + cos(perf_params->th_e * M_PI/180));

    perf_params->C_F_corr = calc_C_F_corr(perf_params->p0, perf_params->p_e,
                                          inputs->p_a, perf_params->Ae_At,
                                          inputs->gamma, perf_params->lm_div);
    perf_params->Ft_corr = perf_params->C_F_corr * flow_props->A[inputs->throat_idx]
                           * perf_params->p0;
    perf_params->Isp_corr = perf_params->C_F_corr * perf_params->cstar / G_0;
}

void Rocket_compute_all(Rocket *self)
{
    Rocket_compute_flow(self);
    Rocket_compute_perf(self);
}

void Rocket_write_summary(Rocket *self, FILE *summary)
{
    PerfParams *perf_params = &self->perf_params;
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
    fprintf(summary, "C_F_corr = %#.6g\n", perf_params->C_F_corr);
    fprintf(summary, "Ft_corr = %#.6g N\n", perf_params->Ft_corr);
    fprintf(summary, "Isp_corr = %#.6g s\n", perf_params->Isp_corr);
    fprintf(summary, "lambda_d = %#.6g\n", perf_params->lm_div);
    fprintf(summary, "theta_e = %#.6g s\n", perf_params->th_e);

    switch (self->flow_props.expansion.flavor) {
        case SUBSONIC:
            fputs("Flow condition = unchoked\n", summary);
            break;
        case SONIC_THROAT:
            fputs("Flow condition = choked subsonic isentropic\n", summary);
            break;
        case NORM_SHK_EXIT:
            fputs("Flow condition = overexpanded w/ normal shock at exit\n", summary);
            fprintf(summary, "Shock x_loc = %.2f\n",
                    self->flow_props.expansion.props.norm_shk_exit.shock_loc);
            break;
        case NORM_SHK_IN:
            fputs("Flow condition = overexpanded w/ normal shock inside\n", summary);
            fprintf(summary, "Shock x_loc = %.2f\n",
                    self->flow_props.expansion.props.norm_shk_in.shock_loc);
            break;
        case OBL_SHK:
            fputs("Flow condition = overexpanded w/ oblique shock\n", summary);
            fprintf(summary, "Shock angle = %#.6g deg\n",
                    self->flow_props.expansion.props.obl_shk.shk_angle);
            break;
        case IDEAL:
            fputs("Flow condition = choked supersonic isentropic\n", summary);
            break;
        case EXP_FAN:
            fputs("Flow condition = underexpanded w/ expansion fan outside\n", summary);
            fprintf(summary, "Turning angle = %#.6g deg\n",
                    self->flow_props.expansion.props.expfan.turn_angle);
            break;
    }
}


int Rocket_write_rawdata(Rocket *self, char *datadir)
{
    int status = 0;
    char *path_buf = calloc(PATH_MAX, sizeof(char));
    RocketInputs *inputs = &self->inputs;
    FlowProperties *flow_props = &self->flow_props;
    size_t n_points = self->inputs.n_points;

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

    free(path_buf);
    return status;
}


int Rocket_make_plots(Rocket *self, char *datadir)
{
    int status = 1;
    char xlabel[] = "normalized x coordinate";
    FILE *gnuplot = popen("gnuplot", "w");
    if (gnuplot == NULL) {
        goto gnuplot_failed;
    }

    gnuplot_plot_array(gnuplot, datadir, "x_position", "mach", xlabel, "Mach number M");
    gnuplot_plot_array(gnuplot, datadir, "x_position", "area", xlabel, "Area A (m^2)");
    gnuplot_plot_array(gnuplot, datadir, "x_position", "flow_velocity", xlabel,
                       "Flow velocity u (m/s)");
    gnuplot_plot_array(gnuplot, datadir, "x_position", "static_density", xlabel,
                       "Static density rho (kg/m^3)");
    gnuplot_plot_with_stagn(gnuplot, datadir, "x_position", "static_pressure",
                            xlabel, "Static pressure p (Pa)",
                            self->perf_params.p0, -5.0, 10.0);
    gnuplot_plot_with_stagn(gnuplot, datadir, "x_position",
                            "static_temperature", xlabel,
                            "Static temperature T (K)",
                            self->perf_params.T0, -5.0, 10.0);
    gnuplot_plot_with_stagn(gnuplot, datadir, "x_position", "specific_enthalpy",
                            xlabel, "Specific enthalpy h (J/kg)",
                            self->perf_params.h0, -5.0, 10.0);
    pclose(gnuplot);
    status = 0;

gnuplot_failed:
    return status;
}
