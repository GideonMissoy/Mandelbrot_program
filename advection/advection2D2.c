#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 31   // number of grid points in x and y directions
#define dt 0.1 // time step
#define Tmax 80.0 // maximum time

int main()
{
    int i, j, n;
    double x, y, t;
    double dx = 1.0, dy = 1.0;
    double ustar = 0.1, k = 0.4, z0 = 0.1;
    double Vx[N][N], Vy[N][N], psi[N][N], z[N];
    double sigma_x = 1.0, sigma_y = 5.0;
    double x0 = 3.0, y0 = 15.0;
    double Q = 100.0/(2*M_PI*sigma_x*sigma_y);
    double u_z, u_z0;

    // Initialize arrays
    for (i = 0; i < N; i++) {
        x = i*dx;
        for (j = 0; j < N; j++) {
            y = j*dy;
            Vx[i][j] = (ustar/k) * log(z[i+1]/z0) - (ustar/k) * log(z[i]/z0); // logarithmic profile
            Vy[i][j] = 0.0;
            psi[i][j] = Q*exp(-pow((x-x0)/sigma_x, 2.0)/2.0 - pow((y-y0)/sigma_y, 2.0)/2.0);
        }
        z[i+1] = (i+1)*dx; // assume constant grid spacing
    }
    z[0] = z[1] - dx; // extrapolate z[0] based on the first two grid points
    u_z0 = (ustar/k) * log(z[1]/z0);

    // Time evolution
    for (n = 0, t = 0.0; t < Tmax; n++, t += dt) {
        // Compute velocities
        for (i = 1; i < N-1; i++) {
            for (j = 1; j < N-1; j++) {
                u_z = (ustar/k) * log(z[i+1]/z[i]); // logarithmic profile
                Vx[i][j] = u_z * psi[i][j+1] - u_z * psi[i][j-1];
                Vy[i][j] = -u_z * psi[i+1][j] + u_z * psi[i-1][j];
            }
        }

        // Update streamfunction
        for (i = 1; i < N-1; i++) {
            for (j = 1; j < N-1; j++) {
                psi[i][j] += dt * (Vx[i][j] - Vy[i][j]);
            }
        }
    }

    // Output results
    FILE *fp;
    fp = fopen("output.dat", "w");
    for (i = 0; i < N; i++) {
        x = i*dx;
        for (j = 0; j < N; j++) {
            y = j*dy;
            fprintf(fp, "%f %f %f\n", x, y, psi[i][j]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);

    return 0;
}

