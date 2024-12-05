#include <stddef.h>
#include <stdio.h>
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
