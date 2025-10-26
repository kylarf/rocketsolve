#include <stddef.h>
#include <stdio.h>
#include <string.h>

#include "fileio.h"

#define BUF_SIZE 65536

size_t count_lines(FILE *file)
{
    char buffer[BUF_SIZE];
    size_t count = 0;
    while (1) {
        size_t bytes_read = fread(buffer, 1, BUF_SIZE, file);
        if (ferror(file)) {
            return -1;
        }
        for (size_t i = 0; i < bytes_read; ++i) {
            if (buffer[i] == '\n') {
                ++count;
            }
        }
        if (feof(file)) {
            break;
        }
    }
    return count;
}

void read_geometry(FILE *geom_file, size_t n_read, double *x_loc, double *nozzle_rad)
{
    // return to position 0 (we've already counted all lines in the file)
    fseek(geom_file, 0, SEEK_SET);
    // skip header line in file
    fscanf(geom_file, "%*[^\n]\n");
    
    for (size_t i = 0; i < n_read; ++i) {
        fscanf(geom_file, "%lf,%lf", &x_loc[i], &nozzle_rad[i]);
    }
}

void read_pressures(FILE *press_file, size_t n_read, double *press_arr)
{
    fseek(press_file, 0, SEEK_SET);

    for (size_t i = 0; i < n_read; ++i) {
        fscanf(press_file, "%lf", &press_arr[i]);
    }
}

int write_array(char *filename, double *arr, size_t N)
{
    int status = 1;

    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        goto fail;
    }

    for (size_t i = 0; i < N; ++i) {
        fprintf(file, "%#.10g\n", arr[i]);
    }

    fclose(file);
    status = 0;
fail:
    return status;
}

int cat_path_write(char *path_buf, char *datadir, char *filename, double *arr, int N)
{
    // make buffer appear to be empty (i.e. ending with null terminator)
    path_buf[0] = '\0';
    // construct path to data file
    strcat(path_buf, datadir);
    //strcat(path_buf, "/");
    strcat(path_buf, filename);
    return write_array(path_buf, arr, N);
}

void gnuplot_plot_array(FILE *gnuplot, char *datadir, char *x_qty_name,
                        char *y_qty_name, char *xlabel, char *ylabel)
{
    fputs("set terminal eps lw 3\n", gnuplot);
    fprintf(gnuplot, "set output '%s%s.eps'\n", datadir, y_qty_name);
    fputs("set size ratio 0.4\n", gnuplot);
    fputs("set yrange [0:*]\n", gnuplot);
    fputs("set offset graph 0, graph 0.1, graph 0.1, graph 0.1\n", gnuplot);
    fprintf(gnuplot, "set ylabel '%s'\n", ylabel);
    fprintf(gnuplot, "set xlabel '%s'\n", xlabel);
    fprintf(gnuplot,
            "plot \"< paste %s%s.dat %s%s.dat\" using 1:2 with lines notitle\n",
            datadir, x_qty_name, datadir, y_qty_name);
}

void gnuplot_plot_with_stagn(FILE *gnuplot, char *datadir, char *x_qty_name,
                             char *y_qty_name, char *xlabel, char *ylabel,
                             double stag, double xmin, double xmax)
{
    fputs("set terminal eps lw 3\n", gnuplot);
    fprintf(gnuplot, "set output '%s%s.eps'\n", datadir, y_qty_name);
    fputs("set size ratio 0.4\n", gnuplot);
    fputs("set yrange [0:*]\n", gnuplot);
    fputs("set offset graph 0, graph 0.1, graph 0.1, graph 0.1\n", gnuplot);
    fprintf(gnuplot, "set ylabel '%s'\n", ylabel);
    fprintf(gnuplot, "set xlabel '%s'\n", xlabel);
    fprintf(gnuplot,
            "plot \"< paste %s%s.dat %s%s.dat\" using 1:2 with lines notitle, "
            "'-' using 1:2 with lines title 'stagnation property'\n",
            datadir, x_qty_name, datadir, y_qty_name);
    fprintf(gnuplot, "%.10g %.10g\n", xmin, stag);
    fprintf(gnuplot, "%.10g %.10g\ne\n", xmax, stag);
}
