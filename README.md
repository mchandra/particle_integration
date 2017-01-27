# FDTD Solver(work in progress)

## Overview

fields.py contains the finite difference EM solver for rectangular grids with uniform cell shapes and interpolator.py contains the function for bilinear interpolation. 

Users can edit **field_diagnostics.py**, **interpolator_diagnostics.py** and **params.py** to run the following tests/tasks:
> **Note:**
> -  Input parameters for the both **field_diagnostics.py** and **interpolator_diagnostics.py** can be edited in **params.py**.


### fields_diagnostics.py (How to use):

The **unit tests** that can be performed with field_diagnostics.py

* Second order error convergence test
* Divergence test

#### Error convergence

*  Edit **N** (test grid size range) near the end of the file as required and run the code.
```python
N = np.array([32, 64, 128, 256])
```

#### Divergence test
movies showing the time evolution of fields and divergence can be made in the following manner.

1. **Comment out** the code for various grid sizes(See near the end of fields_diagnostics.py)
```python
N = np.array([32, 64, 128, 256])
```
and **uncomment the code for fixed grid size** shown below:
```python
#N = np.array([100])
```

2. Dont comment the file writing codes in the script to write data to disk. Make folders named **Ex, Ey, Ez, Bx, By, Bz and div**. The files are saved in the respective folders in h5 format. See code in post.py to read the data. The code below reads data for Ex electric field. 

```python
print('post processing for time_index = ', time_index)
h5f = h5py.File('Ex/solution_'+str(time_index)+'.h5', 'r')
Ex = h5f['Ex/solution_dataset_'+str(time_index)][:]
h5f.close()
```

3. ** Edit post.py** as required to post process the data 

4. Use the code below to make a movie of the images generated. 
```
ffmpeg -f image2 -i point_mass%04d.png -vcodec mpeg4 -mbd rd -trellis 2 -cmp 2 -g 300 -pass 1 -r 25 -b 18000000 movie.mp4
```

### interpolator_diagnostics.py (How to use):

The second order error convergence test can be performed using **interpolator_diagnostics.py**


* Change number of random points taken by editing the line:
```python
  number_random_points = 100
```

*  Change the range for N by editing the lines:
```python
# N = np.array([32, 64, 128, 256, 512, 1024])
N = np.arange(100, 3000, 100)
```

> **Note:**

> -  Use `ghost_cells = 1` for the following code.
> - Read comments before every code segment to see if any thing can be edited as desired.
> - This is a work in progress code and new features and updates will be added on to it.

