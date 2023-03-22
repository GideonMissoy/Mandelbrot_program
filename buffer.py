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
def worker_func(col_buffer, result_buffer):
    while True:
        # Wait for the manager to send a column of the complex plane data
        col_idx = col_buffer.get()
        if col_idx is None:
            break
        
        # Calculate the corresponding column of the Mandelbrot set
        width, height = result_buffer.shape
        for row_idx in range(height):
            c = complex(col_idx / (width/4) - 2, row_idx / (height/4) - 2)
            result_buffer[col_idx, row_idx] = mandelbrot(c)

        # Notify the manager that the results are ready
        result_buffer.flush()
        result_buffer.notify()
    
    # Clean up
    result_buffer.close()

# Define the main function for the manager process
def main():
    # Set up the multiprocessing Queue and shared memory arrays
    col_buffer = mp.Queue()
    result_buffer = np.memmap('/dev/shm/mandelbrot_result.bin', dtype='int32', mode='w+', shape=(1000, 1000))
    
    # Create the worker processes
    num_workers = mp.cpu_count()
    workers = []
    for i in range(num_workers):
        worker = mp.Process(target=worker_func, args=(col_buffer, result_buffer))
        worker.start()
        workers.append(worker)
    
    # Send the columns of the complex plane data to the worker processes
    for col_idx in range(result_buffer.shape[0]):
        col_data = np.random.rand(result_buffer.shape[1]+2) # the buffer size needs to be the size of one column in the complex plane plus two extra values
        col_buffer.put(col_data)
    
    # Tell the worker processes to exit
    for i in range(num_workers):
        col_buffer.put(None

