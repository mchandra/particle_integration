import time
import numpy as np
import arrayfire as af
af.set_backend("opencl")

# This is a benchmark test to compare numpy and arrayfire:

aNumPy = np.random.rand(100, 100)
bNumPy = np.random.rand(100, 100)

np_time_start = time.time()

for i in range(1e6):
  cNumPy = aNumPy + bNumPy

np_time_end     = time.time()
np_time_elapsed = np_time_end - np_time_start

#print("numpy answer is = ", cNumPy)
print("numpy implementation run took time =", np_time_elapsed," seconds")

aArrayFire = af.Array(aNumPy.ctypes.data, aNumPy.shape, aNumPy.dtype.char)
bArrayFire = af.Array(bNumPy.ctypes.data, bNumPy.shape, bNumPy.dtype.char)
cArrayFire = aArrayFire + bArrayFire
cArrayFire.eval()
af.sync()

af_time_start = time.time()

for i in range(1e6):
  cArrayFire = aArrayFire + bArrayFire
  cArrayFire.eval()

af.sync()
af_time_end     = time.time()
af_time_elapsed = af_time_end - af_time_start

#print("arrayfire answer is = ", cArrayFire)
print("arrayfire implementation run took time =", af_time_elapsed," seconds")
