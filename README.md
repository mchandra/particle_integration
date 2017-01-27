# FDTD Solver(work in progress)

## features

fields.py contains the finite difference EM solver for rectangular grids with uniform cell shapes and interpolator.py contains the function for bilinear interpolation. 

Users can edit field_diagnostics.py, interpolator_diagnostics.py and params.py to run the following tests/tasks:

### fields_diagnostics.py (How to use):

For error convergence

1. Edit N (test grid size range) near the end of the file as required and run the code.

For field evolution movies

1. Comment out the code for various grid sizes(See near the end of fields_diagnostics.py) and uncomment the code for fixed grid size.

2. Dont comment the file writing codes in the script to write data to disk. Make folders named Ex, Ey, Ez, Bx, By, Bz and div. The files are saved in the respective folders in h5 format. See code in post.py to read the data. The code below reads data for Ex electric field. 

print('post processing for time_index = ', time_index)
h5f = h5py.File('Ex/solution_'+str(time_index)+'.h5', 'r')
Ex = h5f['Ex/solution_dataset_'+str(time_index)][:]
h5f.close()

3. Edit post.py as required to post process the data 

4. Use the code below to make a movie of the images generated. 

ffmpeg -f image2 -i point_mass%04d.png -vcodec mpeg4 -mbd rd -trellis 2 -cmp 2 -g 300 -pass 1 -r 25 -b 18000000 movie.mp4


### interpolator_diagnostics.py (How to use):

1. Change number of random points taken and test grid ranges as required. 

### Unit tests:

1. Divergence test and error convergence test for fdtd by using fields_diagnostics.py

2. Making movies for field evolution by editing field_diagnostics.py

3. error convergence tests for interpolation using interpolator_diagnostics.py

## How to use the code:

1. Read comments before every code segment to see if any thing can be edited as desired.

2. Download all the files in the current branch and place them in one folder.

3. All the input parameters for the fdtd code can be changed in params.py

4. Use ghostcells = 1 only, for the current solver.

5. This is a work in progress code and new features and updates will be added on to it.
