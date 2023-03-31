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
        z = c;
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

int main(int argc, char *argv[]) {
    int i, j, n, rank, size;
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
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int num_cols_per_process = WIDTH / size;
    int start_col = rank * num_cols_per_process;
    int end_col = start_col + num_cols_per_process - 1;

    double complex *local_plane = malloc((num_cols_per_process + 2) * sizeof(double complex));
    memcpy(local_plane + 1, plane + start_col, num_cols_per_process * sizeof(double complex));

    // Add two extra values to the buffer
    local_plane[0] = plane[start_col-1];
    local_plane[num_cols_per_process+1] = plane[end_col+1];

    int *local_output = malloc(num_cols_per_process * sizeof(int));
    calculate_mandelbrot(start_col, end_col, local_plane, local_output);

    // Send the output of each process to the manager process
    if (rank != 0) {
        MPI_Send(local_output, num_cols_per_process, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } else {
        int *output_ptr_process = output + start_col * HEIGHT;
        memcpy(output_ptr_process, local_output, num_cols_per_process * sizeof(int));

        for (i = 1; i < size; i++) {
            int start_col_i = i * num_cols_per_process;
            int end_col_i = start_col_i + num_cols_per_process - 1;

            int *local_output_i = malloc(num_cols_per_process * sizeof(int));
            MPI_Recv(local_output_i, num_cols_per_process, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            int *output_ptr_process_i = output + start_col_i * HEIGHT;
            memcpy(output_ptr_process_i, local_output_i, num_cols_per_process * sizeof(int));

            free(local_output_i);
        }

        end = clock();

        printf("Execution time: %f seconds\n", (double)(end - start)/CLOCKS_PER_SEC);

        FILE *fp = fopen("mandelbrot.dat", "wb");
	fwrite(output, sizeof(int), WIDTH * HEIGHT, fp);
	fclose(fp);

	free(output);

	return (0);

    }
}
