#!/bin/bash
#SBATCH --job-name=mandelbrot
#SBATCH --output=mandelbrot.out
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --time=00:30:00

# Load any necessary modules
module load mpi

# Run the mandelbrot program using mpirun
mpirun -np 4 mandel1
