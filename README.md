# particle_integration

fields.py contains the finite difference EM solver for rectangular grids with uniform cell shaps and interpolator.py contains the function for bilinear interpolator. Users can edit field_diagnostics.py, interpolator_diagnostics.py and params.py to run the following tests:

1. Divergence test for fdtd.

2. error convergence tests for fdtd and interpolation.



Usage regarding fdtd solver and interpolator:

1. Download all the files in the current branch and place them in one folder.

2. Change all user parameters in params.py

3. All the input parameters for the fdtd code can be changed in params.py

4. Use ghostcells = 1 only, for the current solver.

5. This is a work in progress code and new features and updates will be added on to it.
