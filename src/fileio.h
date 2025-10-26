#ifndef FILEIO_H
#define FILEIO_H

#include "rocket.h"

#include <stdio.h>
#include <stddef.h>

size_t count_lines(FILE *file);

void read_geometry(FILE *geom_file, size_t n_read, double *x_loc, double *nozzle_rad);

void read_pressures(FILE *press_file, size_t n_read, double *press_arr);

int write_array(char *filename, double *arr, size_t N);

int cat_path_write(char *path_buf, char *datadir, char *filename, double *arr, int N);

void gnuplot_plot_array(FILE *gnuplot, char *datadir, char *x_qty_name,
                        char *y_qty_name, char *xlabel, char *ylabel);

void gnuplot_plot_with_stagn(FILE *gnuplot, char *datadir, char *x_qty_name,
                             char *y_qty_name, char *xlabel, char *ylabel,
                             double stag, double xmin, double xmax);

#endif // FILEIO_H
