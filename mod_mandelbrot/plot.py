import matplotlib.pyplot as plt
import numpy as np

# Data for the original program
num_processes = np.array([2, 4, 8, 12, 16, 32])
exec_time_1 = np.array([90.985, 63.48, 34.75, 19.74, 12.31, 8.16])
original_prog = exec_time_1[0] / exec_time_1

# Data for the modified program
num_processes = np.array([2, 4, 8, 12, 16, 32])
exec_time_2 = np.array([61.39, 33.34, 18.41, 10.17, 5.92, 4.61])
mod_prog = exec_time_2[0] / exec_time_2

# Plot the speedup data for the original program
plt.plot(num_processes, original_prog, label='Original')

# Plot the speedup data for the modified program
plt.plot(num_processes, mod_prog, label='Modified')

# Add labels and legend to the plot
plt.title('Speedup vs. Number of MPI Processes')
plt.xlabel('Number of MPI Processes')
plt.ylabel('Parallel Speed-up')
plt.legend()

# Save the plot to a PDF file
plt.savefig('plot.pdf', format='pdf')

