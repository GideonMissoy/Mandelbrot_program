#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

int main(){

	/* Grid properties */
	const int NX=1000;    // Number of x points
  	const int NY=1000;    // Number of y points
  	const float xmin=0.0; // Minimum x value
  	const float xmax=1.0; // Maximum x value
  	const float ymin=0.0; // Minimum y value
  	const float ymax=1.0; // Maximum y value

  	/* Parameters for the Gaussian initial conditions */
	const float x0=0.1;                    // Centre(x)
	const float y0=0.1;                    // Centre(y)
	const float sigmax=0.03;               // Width(x)
	const float sigmay=0.03;               // Width(y)
	const float sigmax2 = sigmax * sigmax; // Width(x) squared
	const float sigmay2 = sigmay * sigmay; // Width(y) squared

	/* Boundary conditions */
  	const float bval_left=0.0;    // Left boudnary value
  	const float bval_right=0.0;   // Right boundary value
  	const float bval_lower=0.0;   // Lower boundary
  	const float bval_upper=0.0;   // Upper bounary
  
  	/* Time stepping parameters */
  	const float CFL=0.9;   // CFL number 
  	const int nsteps=1500; // Number of time steps

  	/* Velocity */
  	const float velx=0.01; // Velocity in x direction
  	const float vely=0.01; // Velocity in y direction

  	/* Arrays to store variables. These have NX+2 elements
        to allow boundary values to be stored at both ends */
  	float x[NX+2];          // x-axis values
  	float y[NX+2];          // y-axis values
  	float u[NX+2][NY+2];    // Array of u values
  	float dudt[NX+2][NY+2]; // Rate of change of u

  	float x2;   // x squared (used to calculate iniital conditions)
  	float y2;   // y squared (used to calculate iniital conditions)
  
	/* Calculate distance between points */
	float dx = (xmax-xmin) / ( (float) NX);
  	float dy = (ymax-ymin) / ( (float) NY);
  
  	/* Calculate time step using the CFL condition */
  	/* The fabs function gives the absolute value in case the velocity is -ve */
  	float dt = CFL / ( (fabs(velx) / dx) + (fabs(vely) / dy) );
  
  	/*** Report information about the calculation ***/
  	printf("Grid spacing dx     = %g\n", dx);
  	printf("Grid spacing dy     = %g\n", dy);
  	printf("CFL number          = %g\n", CFL);
  	printf("Time step           = %g\n", dt);
  	printf("No. of time steps   = %d\n", nsteps);
  	printf("End time            = %g\n", dt*(float) nsteps);
  	printf("Distance advected x = %g\n", velx*dt*(float) nsteps);
  	printf("Distance advected y = %g\n", vely*dt*(float) nsteps);

  	/*** Place x points in parallel for loop ***/
	#pragma omp parallel for
	for (int i = 0; i <= NX + 1; i++)
	{
		x[i] = xmin + ((float) i - 0.5) * dx;
	}

	/* Place y points in parallel for loop */
	#pragma omp parallel for
	for (int j = 0; j <= NY + 1; j++)
	{
		y[j] = ymin + ((float) j - 0.5) * dy;
	}

	/* Set initial conditions */
	#pragma omp parallel for private (x2, y2)
	for (int i = 1; i <= NX; i++)
	{
		x2 = (x[i] - x0) * (x[i] - x0);
		for (int j = 1; j <= NY; j++)
		{
			y2 = (y[j] - y0) * (y[j] - y0);
			u[i][j] = exp(-(x2 / sigmax2) - (y2 / sigmay2));
		}
	}
	/* Apply boundary conditions to initial conditions */
	for (int i = 0; i <= NX + 1; i++)
	{
		u[i][0] = bval_lower;
		u[i][NY + 1] = bval_upper;
	}
	for (int j = 0; j <= NY +1; j++)
	{
		u[0][j] = bval_left;
		u[NX + 1][j] = bval_right;
	}

	/* Time stepping loop */
	for (int n = 1; n <= nsteps; n++)
	{
		/* Calculate rate of change of u using upwinding scheme */
		#pragma omp parallel for
		for (int i = 1; i <= NX; i++)
		{
			for (int j = 1; j <= NY; j++)
			{
				float u_x = (u[i][j] - u[i-1][j]) / dx;
				float u_y = (u[i][j] - u[i][j-1]) / dy;
				if (velx >= 0.0)
				{
					u_x = (u[i+1][j] - u[i][j]) / dx;
				}
				if (vely >= 0.0)
				{
					u_y = (u[i][j+1] - u[i][j]) / dy;
				}
				dudt[i][j] = -velx * (u[i][j] - u[i-1][j]) / dx
					- vely * (u[i][j] - u[i][j-1]) / dy;
			}
		}

		/** Update u **/
		#pragma omp parallel for
		for (int i = 1; i <= NX; i++)
		{
			for (int j = 1; j <= NY; j++)
			{
				u[i][j] = u[i][j] + dudt[i][j] * dt;
			}
		}

		/* Apply boundary conditions */
		#pragma omp parallel for
		for (int i = 0; i <= NX + 1; i++)
		{
			u[i][0] = bval_lower;
			u[i][NY+1] = bval_upper;
		}
		#pragma omp parallel for
		for (int j = 0; j <= NY + 1; j++)
		{
			u[0][j] = bval_left;
			u[NX + 1][j] = bval_right;
		}
	}
	
	/* Output final state of u */
	FILE *gnuplotscript;
	gnuplotscript = fopen("plotscript.gnu", "w");
	fprintf(gnuplotscript, "set pm3d\n");
	fprintf(gnuplotscript, "set view map\n");
	fprintf(gnuplotscript, "set output 'outPut.png'\n");
	fprintf(gnuplotscript, "splot 'final.dat' u 1:2:3\n");
	fclose(gnuplotscript);
	fclose(gp);
	/*call gnuplot to generate a plot */
	system("gnuplot plotscript.gnu");

	return (0);
}
