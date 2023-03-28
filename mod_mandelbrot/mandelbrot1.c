#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <mpi.h>

#define WIDTH 1200
#define HEIGHT 800
#define MAX_ITER 1000

void calculate_mandelbrot(int start_col, int end_col, double complex *plane, int *output) {
    int i, j, k, n;
    double complex c, z;

    for (i = start_col; i <= end_col; i++) {
        c = plane[i];
        z = 0;
        n = 0;
        for (j = 0; j < MAX_ITER; j++) {
            z = z*z + c;
            if (cabs(z) > 2) {
                n = j;
                break;
            }
        }
        output[i - start_col] = n;
    }
}

int main(int argc, char **argv) {
    int i, j, n;
    double complex plane[WIDTH];
    int *output = malloc(WIDTH * HEIGHT * sizeof(int));
    clock_t start, end;

    double x_min = -2.5;
    double x_max = 1;
    double y_min = -1;
    double y_max = 1;

    for (i = 0; i < WIDTH; i++) {
        double x = x_min + (x_max - x_min)/WIDTH * i;
        plane[i] = x + y_min * I;
    }

    MPI_Init(&argc, &argv);

    int num_procs, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_cols_per_proc = WIDTH / num_procs;
    int start_col = rank * num_cols_per_proc;
    int end_col = start_col + num_cols_per_proc - 1;

    double complex *local_plane = malloc((num_cols_per_proc + 2) * sizeof(double complex));
    memcpy(local_plane + 1, plane + start_col, num_cols_per_proc * sizeof(double complex));

    // Add two extra values to the buffer
    local_plane[0] = plane[start_col-1];
    local_plane[num_cols_per_proc+1] = plane[end_col+1];

    int *local_output = malloc(num_cols_per_proc * sizeof(int));

    start = clock();

    calculate_mandelbrot(start_col, end_col, local_plane, local_output);

    MPI_Gather(local_output, num_cols_per_proc, MPI_INT, output, num_cols_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

    end = clock();

    if (rank == 0) {
        printf("Execution time: %f seconds\n", (double)(end - start)/CLOCKS_PER_SEC);

        FILE *fp = fopen("mandelbrot.pgm", "w");
        fprintf(fp, "P2\n%d %d\n%d\n", WIDTH, HEIGHT, MAX_ITER);

        for (j = 0; j < HEIGHT; j++) {
            for (i = 0; i < WIDTH; i++) {
                n = output[i * HEIGHT + j];
                fprintf(fp, "%d ", n);
            }
            fprintf(fp, "\n");
        }

        fclose(fp);
        free(output);
    }

    free(local_plane);
    free(local_output);

    MPI_Finalize();

    return 0;
}
