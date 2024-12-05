#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <sys/syslimits.h>
#include "nozzle.h"
#include "solver.h"
#include "rocket.h"
#include "fileio.h"
#include "stdatm.h"

int main(int argc, char *argv[])
{
    int status = 1;
    FILE *geom_file, *input_file, *summary_file, *pressure_file;

    if (argc < 5) {
        fprintf(stderr, "usage: %s [input-file] [geometry-file] [data-dir] "
                "[pressure-file]\n", argv[0]);
        goto arg_err;
    }

    char *path_buf = calloc(PATH_MAX, sizeof(char));

    pressure_file = fopen(argv[4], "r");
    if (pressure_file == NULL) {
        fprintf(stderr, "Failed to open pressure file: %s\n", argv[4]);
        goto pressure_failed;
    }
    size_t n_press = count_lines(pressure_file);
    double *pressures = calloc(n_press, sizeof(double));
    read_pressures(pressure_file, n_press, pressures);

    // allocate new Rocket for each pressure
    Rocket **rockets_press = malloc(n_press * sizeof(Rocket *));

    geom_file = fopen(argv[2], "r");
    if (geom_file == NULL) {
        fprintf(stderr, "Failed to open geometry file: %s\n", argv[2]);
        goto geom_failed;
    }

    input_file = fopen(argv[1], "r");
    if (input_file == NULL) {
        fprintf(stderr, "Failed to open input file: %s\n", argv[1]);
        goto input_failed;
    }

    Rocket *rocket_default = Rocket_alloc();

    // initialize first ("main") rocket
    if (Rocket_init(rocket_default, input_file, geom_file) != 0) {
        fprintf(stderr, "Failed to parse input file: "
                "at least 1 missing entry or non-numeric leading chars present.\n");
        goto input_parse_failed;
    }

    printf("Read inputs from '%s'.\n", argv[1]);
    printf("Read geometry at %lu points from '%s'.\n\n",
           rocket_default->inputs.n_points, argv[1]);

    Rocket_compute_all(rocket_default);

    for (size_t i = 0; i < n_press; ++i) {
        rockets_press[i] = Rocket_from_inputs(&rocket_default->inputs);
        rockets_press[i]->inputs.p_a = pressures[i];
        Rocket_compute_all(rockets_press[i]);
    }

    snprintf(path_buf, PATH_MAX, "%s/summarized_properties.txt", argv[3]);
    summary_file = fopen(path_buf, "w");
    if (summary_file == NULL) {
        fprintf(stderr, "Failed to open summary file at: %s\n", path_buf);
        goto summary_failed;
    }

    Rocket_write_summary(rocket_default, summary_file);
    fclose(summary_file);
    printf("Wrote summary to '%s'.\n", path_buf);

    // zero out path buffer
    path_buf[0] = '\0';
    snprintf(path_buf, PATH_MAX, "%s/", argv[3]);

    if (Rocket_write_rawdata(rocket_default, path_buf) != 0) {
        fprintf(stderr, "Failed to write one or more raw data files. Check "
                "directory '%s' exists.\n", argv[3]);
    }
    else {
        printf("Wrote raw data to directory '%s'.\n", argv[3]);
    }

    if (Rocket_make_plots(rocket_default, path_buf) != 0) {
        fprintf(stderr, "Failed to write one or more plot files. Check "
                "directory '%s' exists.\n", argv[3]);
    }
    else {
        printf("Wrote plots to directory '%s'.\n", argv[3]);
    }
    puts("");

    for (size_t i = 0; i < n_press; ++i) {
        // zero out path buffer
        path_buf[0] = '\0';
        snprintf(path_buf, PATH_MAX, "%s/%dPa_summarized_properties.txt", argv[3],
                 (int)rockets_press[i]->inputs.p_a);
        summary_file = fopen(path_buf, "w");
        if (summary_file == NULL) {
            fprintf(stderr, "Failed to open summary file with prefix %dPa at: %s\n",
                    (int)rockets_press[i]->inputs.p_a, path_buf);
            goto summary_failed;
        }
        Rocket_write_summary(rockets_press[i], summary_file);
        fclose(summary_file);
        printf("Wrote summary to '%s'.\n", path_buf);
        
        // zero out path buffer again
        path_buf[0] = '\0';
        snprintf(path_buf, PATH_MAX, "%s/%dPa_", argv[3], (int)rockets_press[i]->inputs.p_a);
        if (Rocket_write_rawdata(rockets_press[i], path_buf) != 0) {
            fprintf(stderr, "Failed to write one or more raw data files. Check "
                    "directory '%s' exists.\n", argv[3]);
        }
        else {
            printf("Wrote raw data to directory '%s' with prefix %dPa.\n",
                   argv[3], (int)rockets_press[i]->inputs.p_a);
        }

        if (Rocket_make_plots(rockets_press[i], path_buf) != 0) {
            fprintf(stderr, "Failed to write one or more plot files. Check "
                    "directory '%s' exists.\n", argv[3]);
        }
        else {
            printf("Wrote plots to directory '%s' with prefix %dPa.\n",
                   argv[3], (int)rockets_press[i]->inputs.p_a);
        }
        puts("");
    }

    puts("#### PART 5: ALTITUDE PLOTTING ####");

    int n_alts = 87;
    double altmax = 86;
    double altstep = altmax / (double)n_alts;
    double *alts = malloc(n_alts*sizeof(double));
    double *I_sp = malloc(n_alts*sizeof(double));
    double *C_F = malloc(n_alts*sizeof(double));
    double p, T, rho;
    for (size_t i = 0; i < n_alts; ++i) {
        stdatm(i*altstep, &p, &T, &rho);
        rocket_default->inputs.p_a = p;
        Rocket_compute_all(rocket_default);
        alts[i] = i*altstep;
        I_sp[i] = rocket_default->perf_params.Isp;
        C_F[i] = rocket_default->perf_params.C_F;
    }

    cat_path_write(path_buf, argv[3], "/Isp.dat", I_sp, n_alts);
    cat_path_write(path_buf, argv[3], "/C_F.dat", C_F, n_alts);
    cat_path_write(path_buf, argv[3], "/alts.dat", alts, n_alts);
    path_buf[0] = '\0';
    snprintf(path_buf, PATH_MAX, "%s/", argv[3]);

    FILE *gnuplot = popen("gnuplot", "w");
    gnuplot_plot_array(gnuplot, path_buf, "alts", "Isp", "Altitude (km)", "Isp (s)");
    gnuplot_plot_array(gnuplot, path_buf, "alts", "C_F", "Altitude (km)", "C_F");

    free(alts);
    free(I_sp);
    free(C_F);
    pclose(gnuplot);

    status = 0;

// cleanup

summary_failed:

    Rocket_free(rocket_default);

input_parse_failed:

    fclose(input_file);

input_failed:

    fclose(geom_file);

geom_failed:

    for (size_t i = 0; i < n_press; ++i) {
        Rocket_free(rockets_press[i]);
    }

pressure_failed:

    free(path_buf);

arg_err:

    return status;
}
