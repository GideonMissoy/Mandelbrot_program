#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define WIDTH 1000
#define HEIGHT 1000
#define MAX_ITERATIONS 1000

int mandelbrot(double x, double y) {
    double real = x;
    double imag = y;
    int iterations = 0;
    while (real*real + imag*imag <= 4 && iterations < MAX_ITERATIONS) {
        double real_new = real*real - imag*imag + x;
        double imag_new = 2*real*imag + y;
        real = real_new;
        imag = imag_new;
        iterations++;
    }
    return iterations;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double x_min = -2, x_max = 1, y_min = -1, y_max = 1;
    double x_step = (x_max - x_min) / WIDTH;
    double y_step = (y_max - y_min) / HEIGHT;

    int rows_per_process = HEIGHT / size;
    int start_row = rank * rows_per_process;
    int end_row = (rank + 1) * rows_per_process - 1;

    if (rank == size - 1) {
        end_row = HEIGHT - 1;
    }

    int columns_per_process = WIDTH / size;
    double buffer[columns_per_process + 2];

    for (int row = start_row; row <= end_row; row++) {
        int index = 0;
        for (int col = rank * columns_per_process; col < (rank + 1) * columns_per_process; col++) {
            double x = x_min + col * x_step;
            double y = y_min + row * y_step;
            int iterations = mandelbrot(x, y);
            buffer[index] = x;
            buffer[index + 1] = y;
            buffer[index + 2] = (double)iterations;
            index += 3;
        }
        MPI_Send(buffer, columns_per_process + 2, MPI_DOUBLE, 0, row, MPI_COMM_WORLD);
    }

    if (rank == 0) {
        int row = 0;
        double result[3];
        for (int i = 0; i < HEIGHT; i++) {
            for (int j = 0; j < size; j++) {
                MPI_Recv(buffer, columns_per_process + 2, MPI_DOUBLE, j, row, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                int index = 0;
                for (int col = j * columns_per_process; col < (j + 1) * columns_per_process; col++) {
                    double x = buffer[index];
                    double y = buffer[index + 1];
                    double iterations = buffer[index + 2];
                    result[0] = x;
                    result[1] = y;
                    result[2] = iterations;
                    printf("%f %f %f\n", result[0], result[1], result[2]);
                    index += 3;
                }
                row++;
            }
        }
    }

    MPI_Finalize();

    return 0;
}

