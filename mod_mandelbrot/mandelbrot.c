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
    double complex *plane = malloc(WIDTH * sizeof(double complex));
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

    int num_cols_per_thread = WIDTH / 4;

    start = clock();

    #pragma omp parallel num_threads(4)
    {
        int thread_num = omp_get_thread_num();
        int start_col = thread_num * num_cols_per_thread;
        int end_col = start_col + num_cols_per_thread - 1;

        int *local_output = output + start_col * HEIGHT;

        calculate_mandelbrot(start_col, end_col, plane, local_output);
    }

    end = clock();

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
    free(plane);

    return 0;
}

