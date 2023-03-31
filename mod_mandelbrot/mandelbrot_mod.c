#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <complex.h>
#include <math.h>

#define MAX_ITER 10000
#define MISSING_VALUE -1

int main(int argc, char *argv[]) {
    int rank, size, N;
    double xmin, xmax, ymin, ymax, x, y, xtemp;
    int i, j, iter, *buffer;
    int *recv_counts, *displs, *recv_buffer;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    N = atoi(argv[1]);
    xmin = atof(argv[2]);
    xmax = atof(argv[3]);
    ymin = atof(argv[4]);
    ymax = atof(argv[5]);

    int cols_per_proc = N / size;
    int extra_cols = N % size;

    int local_cols = cols_per_proc + (rank < extra_cols ? 1 : 0);
    int missing_col = (rank >= extra_cols) ? MISSING_VALUE : 0;

    buffer = (int *) malloc((local_cols + 2) * sizeof(int));
    recv_counts = (int *) malloc(size * sizeof(int));
    displs = (int *) malloc(size * sizeof(int));

    for (i = 0; i < local_cols; i++) {
        x = xmin + ((double)(i + (rank * cols_per_proc) + (rank < extra_cols ? rank : extra_cols)) / (double)N) * (xmax - xmin);

        for (j = 0; j < N; j++) {
            y = ymin + ((double)j / (double)N) * (ymax - ymin);

            xtemp = x;
            iter = 0;
            while (x*x + y*y < 4 && iter < MAX_ITER) {
                xtemp = x*x - y*y + x;
                y = 2*x*y + y;
                x = xtemp;
                iter++;
            }

            buffer[i+1] = iter;
        }

        buffer[0] = missing_col;
        buffer[local_cols + 1] = missing_col;
        MPI_Send(buffer, local_cols + 2, MPI_INT, 0, i, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        recv_buffer = (int *) malloc(N * sizeof(int));
        for (i = 0; i < size; i++) {
            recv_counts[i] = cols_per_proc + (i < extra_cols ? 1 : 0);
            displs[i] = i * cols_per_proc + (i < extra_cols ? i : extra_cols);
        }

        for (j = 0; j < N; j++) {
            for (i = 0; i < size; i++) {
                int offset = i * cols_per_proc + (i < extra_cols ? i : extra_cols);
                MPI_Recv(&recv_buffer[offset], cols_per_proc + (i < extra_cols ? 1 : 0), MPI_INT, i, j, MPI_COMM_WORLD, &status);
            }

            for (i = 0; i < N; i++) {
                printf("%d ", recv_buffer[i]);
            }
            printf("\n");
        }

        free(recv_buffer);
    }
    printf("Execution time: %f seconds\n", (double)(end - start)/CLOCKS_PER_SEC);

    FILE *fp = fopen("mandelbrot.dat", "wb");
    fwrite(output, sizeof(int), WIDTH * HEIGHT, fp);
    fclose(fp);

    free(output);
    free(buffer);
    free(recv_counts);
    free(displs);

    MPI_Finalize();
    return 0;
}

