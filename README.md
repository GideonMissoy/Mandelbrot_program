The Mandelbrot set plotted on the complex plane. The horizontal axis is the real axis
and the vertical axis is the imaginary axis.
In the program the number of iterations at each point in the complex plane is stored as a two
dimensional array called nIter. An element of the array is accessed using nIter[i][j] where the
index i refers to the real axis and the index j refers to the imaginary axis. In the manager-worker
version of the program the manger process tells the worker processes which value of the i index to
work on. The worker processes loop over j to calculate a column of values in the complex plane
(this is done in the calc_vals function). When they have completed one column they send a
message to the manager process to request another i value to work on. Workers keep requesting
work until all the values of i have been handed out. Each worker process stores the results it has
calculated, and the results are collated at the end of the calculation using a call to MPI_Reduce
(this is done in the do_communication function).
