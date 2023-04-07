#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>

/*********************************************************************
                      Main function
**********************************************************************/

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

  /*** Place x points in the middle of the cell ***/
  /* LOOP 1 */
/* the loop only has NX+2 iterations, which is a relatively small number, so parallelizing it may not provide any significant performance improvement. */
  for (int i=0; i<NX+2; i++){
    x[i] = ( (float) i - 0.5) * dx;
  }

  /*** Place y points in the middle of the cell ***/
  /* LOOP 2 */
  #pragma omp parallel for
    for (int j=0; j<NY+2; j++){
      y[j] = ( (float) j - 0.5) * dy;
    }


  /*** Set up Gaussian initial conditions ***/
  /* LOOP 3 */
  #pragma omp parallel for
for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
        x2      = (x[i]-x0) * (x[i]-x0);
        y2      = (y[j]-y0) * (y[j]-y0);
        u[i][j] = exp( -1.0 * ( (x2/(2.0*sigmax2)) + (y2/(2.0*sigmay2)) ) );
    }
}


  /*** Write array of initial u values out to file ***/
/* The loop cannot be parallelised, there are no dependencies between the iterations.
the overhead of parallelization may outweigh any potential benefits. */
  FILE *initialfile;
  initialfile = fopen("initial.dat", "w");
  /* LOOP 4 */
  for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
      fprintf(initialfile, "%g %g %g\n", x[i], y[j], u[i][j]);
    }
  }
  fclose(initialfile);
  /*** Update solution by looping over time steps ***/
  /* LOOP 5 */
  for (int m=0; m<nsteps; m++){
    
    /*** Apply boundary conditions at u[0][:] and u[NX+1][:] ***/
    /* LOOP 6 */
/* The loop cannot be parallelised. It is unlikely that parallelizing this loop would provide a significant speedup since there are only NY+2 iteration. */
    for (int j=0; j<NY+2; j++){
      u[0][j]    = bval_left;
      u[NX+1][j] = bval_right;
    }

    /*** Apply boundary conditions at u[:][0] and u[:][NY+1] ***/
    /* LOOP 7 */
/* The loop cannot be parallelised since there are only NX+2 iterations,
there are no data dependencies between the iterations that can be effectively exploited. */
    for (int i=0; i<NX+2; i++){
      u[i][0]    = bval_lower;
      u[i][NY+1] = bval_upper;
    }
    
    /*** Calculate rate of change of u using leftward difference ***/
    /* Loop over points in the domain but not boundary values */
    /* LOOP 8 */
    #pragma omp parallel for
      for (int i=1; i<NX+1; i++){
       for (int j=1; j<NY+1; j++){
         dudt[i][j] = -velx * (u[i][j] - u[i-1][j]) / dx
                 - vely * (u[i][j] - u[i][j-1]) / dy;
       }
     }

    
    /*** Update u from t to t+dt ***/
    /* Loop over points in the domain but not boundary values */
    /* LOOP 9 */
    #pragma omp parallel for
      for (int i=1; i<NX+1; i++){
        for (int j=1; j<NY+1; j++){
          u[i][j] = u[i][j] + dudt[i][j] * dt;
        }
      }

    
  } // time loop
  
  /*** Write array of final u values out to file ***/
    /* Write array of final u values out to file */
  FILE *finalfile;
  finalfile = fopen("final.dat", "w");
  for (int i=0; i<NX+2; i++){
    for (int j=0; j<NY+2; j++){
      fprintf(finalfile, "%g %g %g\n", x[i], y[j], u[i][j]);
    }
    fprintf(finalfile, "\n"); // Blank line to indicate end of row
  }
  fclose(finalfile);

  /* Generate a gnuplot script */
  FILE *gnuplotscript;
  gnuplotscript = fopen("plotscript.gnu", "w");
  fprintf(gnuplotscript, "set pm3d\n");
  fprintf(gnuplotscript, "set view map\n");
  fprintf(gnuplotscript, "splot 'final.dat' u 1:2:3\n");
  fclose(gnuplotscript);

  /* Call gnuplot to generate a plot */
  system("gnuplot plotscript.gnu");


  return 0;
}
