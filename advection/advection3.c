#include <stdio.h>
#include <math.h>

#define Nx 301
#define Ny 301
#define pi M_PI
#define xmax 30.0
#define ymax 30.0
#define dx (xmax/(Nx-1))
#define dy (ymax/(Ny-1))
#define nt 800
#define dt 0.0002
#define alpha 0.1
#define Ustar 0.1
#define kappa 0.4
#define z0 0.01

double f[Nx][Ny], g[Nx][Ny], Vx[Nx][Ny];

int main()
{
    int i, j, n;
    double x, y, t, r, x0, y0, sigmax, sigmay, vx, vy, z, Vh;

    x0 = 3.0;
    y0 = 15.0;
    sigmax = 1.0;
    sigmay = 5.0;

    for(i=0; i<Nx; i++)
    {
        for(j=0; j<Ny; j++)
        {
            x = i*dx;
            y = j*dy;
            r = sqrt((x-x0)*(x-x0)/(sigmax*sigmax) + (y-y0)*(y-y0)/(sigmay*sigmay));
            f[i][j] = exp(-r*r);

            z = y;
            Vh = (Ustar/kappa) * log(z/z0);
            Vx[i][j] = Vh;
        }
    }

    for(n=0; n<nt; n++)
    {
        for(i=1; i<Nx-1; i++)
        {
            for(j=1; j<Ny-1; j++)
            {
                g[i][j] = f[i][j] + alpha*dt/(dx*dx)*(f[i+1][j]-2*f[i][j]+f[i-1][j]) + alpha*dt/(dy*dy)*(f[i][j+1]-2*f[i][j]+f[i][j-1]);
            }
        }
        
        for(i=1; i<Nx-1; i++)
        {
            for(j=1; j<Ny-1; j++)
            {
                f[i][j] = g[i][j] - Vx[i][j]*dt/dx*(g[i][j]-g[i-1][j]) - vy*dt/dy*(g[i][j]-g[i][j-1]);
            }
        }
    }

    FILE *file;
    file = fopen("vert_shear.dat", "w");

    for(i=0; i<Nx; i++)
    {
        for(j=0; j<Ny; j++)
        {
            x = i*dx;
            y = j*dy;
            fprintf(file, "%f %f %f\n", x, y, f[i][j]);
        }
        fprintf(file, "\n");
    }

    fclose(file);

    return 0;
}

