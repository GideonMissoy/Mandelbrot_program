#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <omp.h>

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

int main() {
    int i, j, n;
    double complex plane[WIDTH];
    int output[WIDTH][HEIGHT];
    clock_t start, end;

    double x_min = -2.5;
    double x_max = 1;
    double y_min = -1;
    double y_max = 1;

    for (i = 0; i < WIDTH; i++) {
        double x = x_min + (x_max - x_min)/WIDTH * i;
        plane[i] = x + y_min * I;
    }

    int *output_ptr = &output[0][0];
    int num_cols_per_thread = WIDTH / 4;

    start = clock();

    #pragma omp parallel num_threads(4)
    {
        int thread_num = omp_get_thread_num();
        int start_col = thread_num * num_cols_per_thread;
        int end_col = start_col + num_cols_per_thread - 1;

        double complex local_plane[num_cols_per_thread+2];
        memcpy(local_plane + 1, plane + start_col, num_cols_per_thread * sizeof(double complex));

        // Add two extra values to the buffer
        local_plane[0] = plane[start_col-1];
        local_plane[num_cols_per_thread+1] = plane[end_col+1];

        int local_output[num_cols_per_thread];
        calculate_mandelbrot(start_col, end_col, local_plane, local_output);

        // Copy the output of each thread back to the main output array
        int *output_ptr_thread = output_ptr + start_col * HEIGHT;
        memcpy(output_ptr_thread, local_output, num_cols_per_thread * sizeof(int));
    }

    end = clock();

    printf("Execution time: %f seconds\n", (double)(end - start)/CLOCKS_PER_SEC);

    FILE *fp = fopen("mandelbrot.pgm", "w");
    fprintf(fp, "P2\n%d %d\n%d\n", WIDTH, HEIGHT, MAX_ITER);

    for (j = 0; j < HEIGHT; j++) {
        for (i = 0; i < WIDTH; i++) {
            n = output[i][j];
            fprintf(fp, "%d ", n);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);

    return 0;
}
