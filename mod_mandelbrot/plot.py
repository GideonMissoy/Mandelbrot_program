import matplotlib.pyplot as plt
import numpy as np

# Enter the time taken by original program for 1, 2, 4, 8, 16 and 32 processes
orig_times = [299.38, 170.07, 93.83, 52.23, 32.10, 23.97]

# Enter the time taken by modified program for 1, 2, 4, 8, 16 and 32 processes
mod_times = [335.47, 182.56, 99.87, 55.45, 32.50, 24.71]

# Calculate speed up for modified program
speed_up = [orig_times[0]/time for time in mod_times]

# Plot the graph
x = np.array([1, 2, 4, 8, 16, 32])
plt.plot(x, x, 'r--', label='Ideal speed up')
plt.plot(x, x/speed_up, '-o', label='Modified program')
plt.plot(x, x/orig_times, '-o', label='Original program')
plt.xlabel('Number of MPI processes')
plt.ylabel('Speed up')
plt.legend(loc='best')
plt.show()

