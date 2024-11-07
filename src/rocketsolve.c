#include <stdio.h>
#include <string.h>
#include "solver.h"
#include "rocket.h"

int main(int argc, char *argv[])
{
    int status = 1;
    FILE *geom_file, *input_file, *summary_file;

    if (argc < 4) {
        fprintf(stderr, "usage: %s [input-file] [geometry-file] [data-dir]\n", argv[0]);
        goto arg_err;
    }

    char *path_buf = calloc(4096, sizeof(char));
    Rocket *rocket = Rocket_alloc();

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

    if (Rocket_init(rocket, input_file, geom_file) != 0) {
        fprintf(stderr, "Failed to parse input file: "
                "at least 1 missing entry or non-numeric leading chars present.\n");
        goto input_parse_failed;
    }

    printf("Read inputs from '%s'.\n", argv[1]);
    printf("Read geometry at %lu points from '%s'.\n", rocket->inputs->n_points, argv[1]);

    Rocket_compute_flow(rocket);

    Rocket_compute_perf(rocket);

    strcat(path_buf, argv[3]);
    strcat(path_buf, "/summarized_properties.txt");
    summary_file = fopen(path_buf, "w");
    if (summary_file == NULL) {
        fprintf(stderr, "Failed to open summary file at: %s\n", path_buf);
        goto summary_failed;
    }

    Rocket_write_summary(rocket, summary_file);
    printf("Wrote summary to '%s'.\n", path_buf);

    if (Rocket_write_rawdata(rocket, argv[3]) != 0) {
        fprintf(stderr, "Failed to write one or more raw data files. Check "
                "directory '%s' exists.\n", argv[3]);
    }
    else {
        printf("Wrote raw data to directory '%s'.\n", argv[3]);
    }

    if (Rocket_make_plots(rocket, argv[3]) != 0) {
        fprintf(stderr, "Failed to write one or more plot files. Check "
                "directory '%s' exists.\n", argv[3]);
    }
    else {
        printf("Wrote plots to directory '%s'.\n", argv[3]);
    }
  
    status = 0;

// cleanup

    fclose(summary_file);

summary_failed:
input_parse_failed:

    fclose(input_file);

input_failed:

    fclose(geom_file);

geom_failed:

    Rocket_free(rocket);
    free(path_buf);

arg_err:

    return status;
}
