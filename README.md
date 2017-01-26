# FDTD Solver(work in progress)

fields.py contains the finite difference EM solver for rectangular grids with uniform cell shaps and interpolator.py contains the function for bilinear interpolator. Users can edit field_diagnostics.py, interpolator_diagnostics.py and params.py to run the following tests/tasks:

1. Divergence test and error convergence test for fdtd by using fields_diagnostics.py

2. Making movies for field evoulution by editing field_diagnostics.py

3. error convergence tests for interpolation using 



Usage regarding fdtd solver and interpolator:

1. Read comments before every code segment to see if any thing can be edited as desired.

2. Download all the files in the current branch and place them in one folder.

3. All the input parameters for the fdtd code can be changed in params.py

4. Use ghostcells = 1 only, for the current solver.

5. This is a work in progress code and new features and updates will be added on to it.
