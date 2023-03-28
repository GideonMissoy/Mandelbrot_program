import matplotlib.pyplot as plt
import numpy as np

num_processes = np.array([2, 4, 8, 12, 16, 32])
exec_time = np.array([61.39, 33.34, 18.41, 10.17, 5.92, 4.61])
speedup = exec_time[0] / exec_time

plt.plot(num_processes, speedup)
plt.title('Speedup vs. Number of MPI Processes')
plt.xlabel('Number of MPI Processes')
plt.ylabel('Speedup')

plt.savefig('plot.pdf', format='pdf')
