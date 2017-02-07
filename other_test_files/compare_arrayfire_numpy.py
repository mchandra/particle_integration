import time
import numpy as np
import arrayfire as af
af.set_backend("opencl")

# This is a benchmark test to compare numpy and arrayfire:

np_time_start = time.time()

a = np.array([1, 2, 3, 4, 5])
b = np.array([6, 7, 8, 9, 10])
c = a + b

np_time_end     = time.time()
np_time_elapsed = np_time_end - np_time_start

print("numpy answer is = ", c)
print("numpy implementation run took time =", np_time_elapsed," seconds")

af_time_start = time.time()

a = af.Array([1, 2, 3, 4, 5])
b = af.Array([6, 7, 8, 9, 10])
c = a + b

af_time_end     = time.time()
af_time_elapsed = af_time_end - af_time_start

print("arrayfire answer is = ", c)
print("arrayfire implementation run took time =", af_time_elapsed," seconds")
