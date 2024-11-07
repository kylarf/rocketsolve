#ifndef FILEIO_H
#define FILEIO_H

#include <stdio.h>
#include "rocket.h"

size_t count_lines(FILE *file);

void read_geometry(FILE *geom_file, size_t n_read, double *x_loc, double *nozzle_rad);

int write_array(char *filename, double *arr, size_t N);

#endif // FILEIO_H
