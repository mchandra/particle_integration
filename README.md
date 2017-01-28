# FDTD Solver(work in progress)

## Overview

fields.py contains the finite difference EM solver for rectangular grids with uniform cell shapes and interpolator.py contains the function for bilinear interpolation.

Users can edit **field_diagnostics.py**, **interpolator_diagnostics.py** and **params.py** to run the following tests/tasks:
> **Note:**
> -  Input parameters for the both **field_diagnostics.py** and **interpolator_diagnostics.py**  collectively can be edited in **params.py**.


### fields_diagnostics.py (How to use):

The **unit tests** that can be performed with field_diagnostics.py

* Second order error convergence test
* Divergence test

#### Error convergence

*  Edit **N** (test grid size range) found near the end of the file as required and run the code.
```python
N = np.array([32, 64, 128, 256])
```

#### Divergence test
movies showing the time evolution of fields and divergence can be made in the following manner.

* **Comment out** the code for various grid sizes(See near the end of fields_diagnostics.py)
```python
N = np.array([32, 64, 128, 256])
```
and **uncomment the code for fixed grid size** shown below:
```python
#N = np.array([100])
```

* Uncomment the following set of lines used for writing data to disk. They have been commented out by default.

```
    h5f = h5py.File('Ex/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('Ex/solution_dataset_'+str(time_index), data=Ex)
    h5f.close()

    h5f = h5py.File('Ey/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('Ey/solution_dataset_'+str(time_index), data=Ey)
    h5f.close()

    h5f = h5py.File('Ez/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('Ez/solution_dataset_'+str(time_index), data=Ez)
    h5f.close()

    h5f = h5py.File('Bx/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('Bx/solution_dataset_'+str(time_index), data=Bx)
    h5f.close()

    h5f = h5py.File('By/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('By/solution_dataset_'+str(time_index), data=By)
    h5f.close()

    h5f = h5py.File('Bz/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('Bz/solution_dataset_'+str(time_index), data=Bz)
    h5f.close()

    h5f = h5py.File('div/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('div/solution_dataset_'+str(time_index), data=div_B)
    h5f.close()
```

Make sure there are folders named **Ex, Ey, Ez, Bx, By, Bz and div** in your working directory. The data generated from the code will be saved in `.h5` format. See post.py to see scripts for reading the data. The code below shows a sample script for reading the data generated for `Ex` electric field.

```python
print('post processing for time_index = ', time_index)
h5f = h5py.File('Ex/solution_'+str(time_index)+'.h5', 'r')
Ex = h5f['Ex/solution_dataset_'+str(time_index)][:]
h5f.close()
```

* ** Edit post.py** as required to post process the data

* Use the code below to make a movie of the images generated.

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
*  Edit  the following set of statements to plot the desired variables.

```
pl.loglog(N, ErrorNEz, '-o', label='$E_z$ ')
pl.legend()
pl.loglog(N, ErrorNBx, '-o', label='$B_x$ ')
pl.legend()
pl.loglog(N, ErrorNBy, '-o', label='$B_y$ ')
pl.legend()
pl.loglog(N, ErrorNBz, label='$B_z$ ')
pl.legend()
pl.loglog(N, ErrorNEx, label='$E_x$ ')
pl.legend()
pl.loglog(N, ErrorNEy, label='$E_y$ ')
pl.legend()
pl.loglog(N, 1.5 * (N ** -1.999), '--', color='black', lw=2, label=' $O(N^{-2})$ ')
pl.legend()
pl.title('$\mathrm{Convergence\; plot}$ ')
pl.xlabel('$\mathrm{N}$')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
pl.show()
pl.clf()
```

> **Important Note:**

> -  Use `ghost_cells = 1` for all the files.
> - Read comments before every code segment to see if any thing can be edited as desired.
> - This is a work in progress code and new features and updates will be added on to it.
