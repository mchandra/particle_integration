import numpy as np

""" Equations for mode 1 fdtd"""

# dEz/dt = dBy/dx - dBx/dy
# dBx/dt = -dEz/dy
# dBy/dt = +dEz/dx
# div_B = dBx/dx + dBy/dy


""" Equations for mode 2 fdtd"""

# dBz/dt = - ( dEy/dx - dEx/dy )
# dEx/dt = + dBz/dy
# dEy/dt = - dBz/dx
# div_B = dBz/dz

""" Alignment of the spatial grids (Convention chosen)"""

# Alignment of the x direction

"""
mode 1:

Ez (x_center)
------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o---
            |
            |
---------starts here--------

Bx (x_center)
------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o---
            |
            |
---------starts here--------

By (x_right)
----------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------
                |
                |
-----------starts here--------

mode :

Bz (x_right)
----------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------
                |
                |
-----------starts here--------

Ex (x_right)
----------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------
                |
                |
-----------starts here--------

Ey (x_center)
------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o---
            |
            |
---------starts here--------

"""

# Alignment of the y direction

"""
mode 1:

Ez (y_center)
------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o---
            |
            |
---------starts here--------

Bx (y_right)
----------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------
                |
                |
-----------starts here--------

By (y_center)
------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o---
            |
            |
---------starts here--------

mode :

Bz (y_right)
----------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------
                |
                |
-----------starts here--------

Ex (y_center)
------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o---
            |
            |
---------starts here--------
Ey (y_right)
----------------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------o-------
                |
                |
-----------starts here--------
"""


"""  Setting number of ghost cells  """


def fdtd(Ex, Ey, Ez, Bx, By, Bz, c, Lx, Ly, ghost_cells, Jx, Jy, Jz):


  Ez_updated, Bx_updated, By_updated = mode1_fdtd(Ez, Bx, By, Lx, Ly, c, ghost_cells, Jx, Jy, Jz)

  Bz_updated, Ex_updated, Ey_updated = mode2_fdtd(Bz, Ex, Ey, Lx, Ly, c, ghost_cells, Jx, Jy, Jz )

  return Ex_updated, Ey_updated, Ez_updated, Bx_updated, By_updated, Bz_updated


""" Equations for mode 1 fdtd (variation along x and y)"""

# dEz/dt = dBy/dx - dBx/dy
# dBx/dt = -dEz/dy
# dBy/dt = +dEz/dx
# div_B = dBx/dx + dBy/dy


def mode1_fdtd( Ez, Bx, By, Lx, Ly, c, ghost_cells, Jx, Jy, Jz ):

  """ Number of divisions in the physical domain"""

  x_number_of_points,  y_number_of_points = Ez.shape

  Nx = x_number_of_points - 2*ghost_cells - 1
  Ny = y_number_of_points - 2*ghost_cells - 1

  Ez_in_function = Ez
  Bx_in_function = Bx
  By_in_function = By

  """Enforcing BC's"""

  Ez_in_function[0, :]                      = Ez_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Ez_in_function[:, 0]                      = Ez_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Ez_in_function[y_number_of_points - 1, :] = Ez_in_function[ghost_cells, :].copy()
  Ez_in_function[:, x_number_of_points - 1] = Ez_in_function[:, ghost_cells].copy()

  Bx_in_function[0, :]                      = Bx_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Bx_in_function[:, 0]                      = Bx_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Bx_in_function[y_number_of_points - 1, :] = Bx_in_function[ghost_cells, :].copy()
  Bx_in_function[:, x_number_of_points - 1] = Bx_in_function[:, ghost_cells].copy()

  By_in_function[0, :]                      = By_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  By_in_function[:, 0]                      = By_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  By_in_function[y_number_of_points - 1, :] = By_in_function[ghost_cells, :].copy()
  By_in_function[:, x_number_of_points - 1] = By_in_function[:, ghost_cells].copy()

  """ Setting division size and time steps"""

  dx = np.float(1 / (Nx))
  dy = np.float(1 / (Ny))
  dt = np.float(dx / (2 * c))


  """ defining variable for convenience """

  dt_by_dx = dt / (dx)
  dt_by_dy = dt / (dy)

  """  Defining index grid for updating the fields  """

  I, J = np.meshgrid(range(ghost_cells, x_number_of_points-ghost_cells), \
                      range(ghost_cells, y_number_of_points-ghost_cells)
                      )

  """  Updating the Electric field  """

  Ez_in_function[I, J] = Ez_in_function[I, J] + ( (dt_by_dx * (By_in_function[I, J] - By_in_function[I, J - 1]))\
                                                 - (dt_by_dy * (Bx_in_function[I, J] - Bx_in_function[I - 1, J]))\
                                                 )

  # dEz/dt = dBy/dx - dBx/dy

  """  Implementing periodic boundary conditions using ghost cells  """

  Ez_in_function[0, :]                      = Ez_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Ez_in_function[:, 0]                      = Ez_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Ez_in_function[y_number_of_points - 1, :] = Ez_in_function[ghost_cells, :].copy()
  Ez_in_function[:, x_number_of_points - 1] = Ez_in_function[:, ghost_cells].copy()

  """  Updating the Magnetic fields   """

  Bx_in_function[I, J] = Bx_in_function[I, J] - (dt_by_dy * (Ez_in_function[I + 1, J] - Ez_in_function[I, J]))

  # dBx/dt = -dEz/dy

  By_in_function[I, J] = By_in_function[I, J] + (dt_by_dx * (Ez_in_function[I, J + 1] - Ez_in_function[I, J]))

  # dBy/dt = +dEz/dx

  """  Implementing periodic boundary conditions using ghost cells  """

  Bx_in_function[0, :]                      = Bx_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Bx_in_function[:, 0]                      = Bx_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Bx_in_function[y_number_of_points - 1, :] = Bx_in_function[ghost_cells, :].copy()
  Bx_in_function[:, x_number_of_points - 1] = Bx_in_function[:, ghost_cells].copy()



  By_in_function[0, :]                      = By_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  By_in_function[:, 0]                      = By_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  By_in_function[y_number_of_points - 1, :] = By_in_function[ghost_cells, :].copy()
  By_in_function[:, x_number_of_points - 1] = By_in_function[:, ghost_cells].copy()

  return Ez_in_function, Bx_in_function, By_in_function



"""-------------------------------------------------End--of--Mode--1-------------------------------------------------"""


"""-------------------------------------------------Start--of--Mode-2------------------------------------------------"""


""" Equations for mode 2 fdtd (variation along x and y)"""


# dBz/dt = - ( dEy/dx - dEx/dy )
# dEx/dt = + dBz/dy
# dEy/dt = - dBz/dx
# div_B = dBz/dz


def mode2_fdtd( Bz, Ex, Ey, Lx, Ly, c, ghost_cells, Jx, Jy, Jz ):

  """ Number of divisions in the physical domain"""

  x_number_of_points,  y_number_of_points = Bz.shape

  Nx = x_number_of_points - 2*ghost_cells-1
  Ny = y_number_of_points - 2*ghost_cells-1

  Bz_in_function = Bz
  Ex_in_function = Ex
  Ey_in_function = Ey

  """Enforcing BC's"""

  Bz_in_function[0, :]                      = Bz_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Bz_in_function[:, 0]                      = Bz_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Bz_in_function[y_number_of_points - 1, :] = Bz_in_function[ghost_cells, :].copy()
  Bz_in_function[:, x_number_of_points - 1] = Bz_in_function[:, ghost_cells].copy()

  Ex_in_function[0, :]                      = Ex_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Ex_in_function[:, 0]                      = Ex_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Ex_in_function[y_number_of_points - 1, :] = Ex_in_function[ghost_cells, :].copy()
  Ex_in_function[:, x_number_of_points - 1] = Ex_in_function[:, ghost_cells].copy()

  Ey_in_function[0, :]                      = Ey_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Ey_in_function[:, 0]                      = Ey_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Ey_in_function[y_number_of_points - 1, :] = Ey_in_function[ghost_cells, :].copy()
  Ey_in_function[:, x_number_of_points - 1] = Ey_in_function[:, ghost_cells].copy()


  """ Setting division size and time steps"""

  dx = np.float(1 / (Nx))
  dy = np.float(1 / (Ny))
  dt = np.float(dx / (2 * c))

  """ defining variable for convenience """

  dt_by_dx = dt / (dx)
  dt_by_dy = dt / (dy)

  """  Defining index grid for updating the fields  """

  I, J = np.meshgrid(range(ghost_cells, x_number_of_points-ghost_cells), \
                      range(ghost_cells, y_number_of_points-ghost_cells)
                      )

  """  Updating the Electric field  """

  Bz_in_function[I, J] = Bz_in_function[I, J] - ((dt_by_dx * (Ey_in_function[I, J] - Ey_in_function[I, J - 1]))\
                                                  - (dt_by_dy * (Ex_in_function[I, J] - Ex_in_function[I - 1, J]))\
                                                 )

  # dBz/dt = - ( dEy/dx - dEx/dy )

  """  Implementing periodic boundary conditions using ghost cells  """

  Bz_in_function[0, :]                      = Bz_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Bz_in_function[:, 0]                      = Bz_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Bz_in_function[y_number_of_points - 1, :] = Bz_in_function[ghost_cells, :].copy()
  Bz_in_function[:, x_number_of_points - 1] = Bz_in_function[:, ghost_cells].copy()

  """  Updating the Magnetic fields   """

  Ex_in_function[I, J] = Ex_in_function[I, J] + (dt_by_dy * (Bz_in_function[I + 1, J] - Bz_in_function[I, J]))

  # dEx/dt = + dBz/dy

  Ey_in_function[I, J] = Ey_in_function[I, J] - (dt_by_dx * (Bz_in_function[I, J + 1] - Bz_in_function[I, J]))

  # dEy/dt = - dBz/dx

  """  Implementing periodic boundary conditions using ghost cells  """

  Ex_in_function[0, :]                      = Ex_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Ex_in_function[:, 0]                      = Ex_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Ex_in_function[y_number_of_points - 1, :] = Ex_in_function[ghost_cells, :].copy()
  Ex_in_function[:, x_number_of_points - 1] = Ex_in_function[:, ghost_cells].copy()

  Ey_in_function[0, :]                      = Ey_in_function[y_number_of_points - 1 - ghost_cells, :].copy()
  Ey_in_function[:, 0]                      = Ey_in_function[:, x_number_of_points - 1 - ghost_cells].copy()
  Ey_in_function[y_number_of_points - 1, :] = Ey_in_function[ghost_cells, :].copy()
  Ey_in_function[:, x_number_of_points - 1] = Ey_in_function[:, ghost_cells].copy()

  return Bz_in_function, Ex_in_function, Ey_in_function

"""-------------------------------------------------End--of--Mode--2-------------------------------------------------"""
