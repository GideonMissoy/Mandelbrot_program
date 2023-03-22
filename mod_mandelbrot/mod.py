import numpy as np
import multiprocessing as mp

# Define the function to calculate the Mandelbrot set
def mandelbrot(c, max_iters=100):
    z = 0
    for i in range(max_iters):
        z = z**2 + c
        if abs(z) > 2:
            return i
    return max_iters

# Define the function for the worker process
def worker_func(col_idx, col_data, result_queue):
    # Calculate the corresponding column of the Mandelbrot set
    width, height = col_data.shape
    result_col = np.zeros(height, dtype='int32')
    for row_idx in range(height):
        c = complex(col_idx / (width/4) - 2, row_idx / (height/4) - 2)
        result_col[row_idx] = mandelbrot(c)

    # Send the results back to the manager process
    result_queue.put((col_idx, result_col))

# Define the main function for the manager process
def main():
    # Set up the multiprocessing Queue and shared memory arrays
    result_buffer = np.memmap('/dev/shm/mandelbrot_result.bin', dtype='int32', mode='w+', shape=(1000, 1000))

    # Create the worker processes
    num_workers = mp.cpu_count()
    workers = []
    result_queue = mp.Queue()
    for i in range(num_workers):
        worker = mp.Process(target=worker_func, args=(i, result_buffer[i:i+1,:], result_queue))
        worker.start()
        workers.append(worker)

    # Receive the results from the worker processes
    for i in range(result_buffer.shape[0]):
        col_idx, result_col = result_queue.get()
        result_buffer[col_idx, :] = result_col

    # Tell the worker processes to exit
    for i in range(num_workers):
        workers[i].join()

if __name__ == '__main__':
    main()

