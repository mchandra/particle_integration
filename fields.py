import numpy as np

""" Equations for mode 1 fdtd"""

# dEz/dt = dBy/dx - dBx/dy
# dBx/dt = -dEz/dy
# dBy/dt = +dEz/dx
# div_B  = dBx/dx + dBy/dy


""" Equations for mode 2 fdtd"""

# dBz/dt = - ( dEy/dx - dEx/dy )
# dEx/dt = + dBz/dy
# dEy/dt = - dBz/dx
# div_B  = dBz/dz

"""
Notes for periodic boundary conditions:
for [0, Lx] domain use periodic BC's such that last point in the physical domain coincides with the first point
for [0, Lx) domain use periodic BC's such that the ghost point after the last physical point coincides with the first physical point
"""


""" Alignment of the spatial grids for the fields(Convention chosen)"""

"""


positive y axis going down due to matrix representation in computers
positive x axis -------------> going right

Let the domain be [0,1]
Sample grid with one ghost cell at each end and the physical domain containing only 2 points
Here dx = 1, dx/2 = 0.5

Let the grids for the example case be denoted be:

x_center = [-1, 0, 1, 2]
y_center = [-1, 0, 1, 2]

x_center[0] and x_center[3] are the ghost points and x_center[1] and x_center[2] are the physical points
y_center[0] and y_center[3] are the ghost points and y_center[1] and y_center[2] are the physical points


x_right  = [-0.5, 0.5, 1.5, 2.5]
y_top    = [-0.5, 0.5, 1.5, 2.5]

x_right[0] and x_right[3] are the ghost points and x_right[1] and x_right[2] are the physical points
y_top[0] and y_top[3] are the ghost points and y_top[1] and y_top[2] are the physical points

This can be seen visually with the below presented schematic



+++++++++p11--------p12--------p13--------p14--------p15--------p16--------p17--------p18+++++++++++++++++++++++++++++++++++++++++++++

p11 = (x_center[0], y_center[0]) (ghost point), p13 = (x_center[1], y_center[0]) (physical point)
p15 = (x_center[2], y_center[0]) (physical point),p17 = (x_center[3], y_center[0]) (ghost point)
p12 = (x_right[0], y_center[0]) (ghost point), p14 = (x_right[1], y_center[0]) (physical point)
p16 = (x_right[2], y_center[0]) (physical point), p18 = (x_right[3], y_center[0]) (ghost point)

+++++++++p21--------p22--------p23--------p24--------p25--------p26--------p27--------p28+++++++++++++++++++++++++++++++++++++++++++++

p21 = (x_center[0], y_top[0]) (ghost point), p23 = (x_center[1], y_top[0]) (physical point)
p25 = (x_center[2], y_top[0]) (physical point), p27 = (x_center[3], y_top[0]) (ghost point)
p22 = (x_right[0], y_top[0]) (ghost point), p24 = (x_right[1], y_top[0]) (physical point)
p26 = (x_right[2], y_top[0]) (physical point), p28 = (x_right[3], y_top[0]) (ghost point)

+++++++++p31--------p32--------p33--------p34--------p35--------p36--------p37--------p38+++++++++++++++++++++++++++++++

p31 = (x_center[0], y_center[1]) (ghost point), p33 = (x_center[1], y_center[1]) (physical point)
p35 = (x_center[2], y_center[1]) (physical point), p37 = (x_center[3], y_center[1]) (ghost point)
p32 = (x_right[0], y_center[1]) (ghost point), p34 = (x_right[1], y_center[1]) (physical point)
p36 = (x_right[2], y_center[1]) (physical point), p38 = (x_right[3], y_center[1]) (ghost point)

+++++++++p41--------p42--------p43--------p44--------p45--------p46--------p47--------p48+++++++++++++++++++++++++++++++

p41 = (x_center[0], y_top[1]) (ghost point), p43 = (x_center[1], y_top[1]) (physical point)
p45 = (x_center[2], y_top[1]) (physical point), p47 = (x_center[3], y_top[1]) (ghost point)
p42 = (x_right[0], y_top[1]) (ghost point), p44 = (x_right[1], y_top[1]) (physical point)
p46 = (x_right[2], y_top[1]) (physical point), p48 = (x_right[3], y_top[1]) (ghost point)

+++++++++p51--------p52--------p53--------p54--------p55--------p56--------p57--------p58+++++++++++++++++++++++++++++++

And So on ................




+++++++++p61--------p62--------p63--------p64--------p65--------p66--------p67--------p68+++++++++++++++++++++++++++++++

And So on ................




+++++++++p71--------p72--------p73--------p74--------p75--------p76--------p77--------p78+++++++++++++++++++++++++++++++

And So on ................




+++++++++p71--------p72--------p73--------p74--------p75--------p76--------p77--------p88+++++++++++++++++++++++++++++++


Now the fields aligned in x and y direction along with the following grids:

Ex = (x_right, y_center )
Ey = (x_center, y_top   )
Ez = (x_center, y_center)
Bx = (x_center, y_top   )
By = (x_right, y_center )
Bz = (x_right, y_top    )




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

""" Function for enforcing periodic boundary conditions on a field"""

def periodic(field, x_points, y_points, ghost_cells):

    field[0, :]            = field[y_points - 2 - ghost_cells, :].copy()
    field[:, 0]            = field[:, x_points - 2 - ghost_cells].copy()
    field[y_points - 1, :] = field[ghost_cells + 1, :].copy()
    field[:, x_points - 1] = field[:, ghost_cells + 1].copy()

    return field



def mode1_fdtd( Ez, Bx, By, Lx, Ly, c, ghost_cells, Jx, Jy, Jz ):

  """ Number of divisions in the physical domain"""

  x_number_of_points,  y_number_of_points = Ez.shape

  Nx = x_number_of_points - 2*ghost_cells - 1
  Ny = y_number_of_points - 2*ghost_cells - 1

  Ez_in_function = Ez
  Bx_in_function = Bx
  By_in_function = By

  """Enforcing BC's"""

  Ez_in_function = periodic(Ez_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  Bx_in_function = periodic(Bx_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  By_in_function = periodic(By_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  """ Setting division size and time steps"""

  dx = np.float(Lx / (Nx))
  dy = np.float(Ly / (Ny))
  dt = np.float(dx / (2 * c))


  """ defining variable for convenience """

  dt_by_dx = dt / (dx)
  dt_by_dy = dt / (dy)

  """  Defining index grid for updating the fields  """

  I, J = np.meshgrid(range(ghost_cells, x_number_of_points-ghost_cells), \
                      range(ghost_cells, y_number_of_points-ghost_cells)\
                    )

  """  Updating the Electric field  """

  Ez_in_function[I, J] = Ez_in_function[I, J] + ( (dt_by_dx * (By_in_function[I, J] - By_in_function[I, J - 1]))\
                                                 - (dt_by_dy * (Bx_in_function[I, J] - Bx_in_function[I - 1, J]))\
                                                )

  # dEz/dt = dBy/dx - dBx/dy

  """  Implementing periodic boundary conditions using ghost cells  """

  Ez_in_function = periodic(Ez_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  """  Updating the Magnetic fields   """

  Bx_in_function[I, J] = Bx_in_function[I, J] - (dt_by_dy * (Ez_in_function[I + 1, J] - Ez_in_function[I, J]))

  # dBx/dt = -dEz/dy

  By_in_function[I, J] = By_in_function[I, J] + (dt_by_dx * (Ez_in_function[I, J + 1] - Ez_in_function[I, J]))

  # dBy/dt = +dEz/dx

  """  Implementing periodic boundary conditions using ghost cells  """

  Bx_in_function = periodic(Bx_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  By_in_function = periodic(By_in_function, x_number_of_points, y_number_of_points, ghost_cells)

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

  """Enforcing periodic BC's"""

  Bz_in_function = periodic(Bz_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  Ex_in_function = periodic(Ex_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  Ey_in_function = periodic(Ey_in_function, x_number_of_points, y_number_of_points, ghost_cells)


  """ Setting division size and time steps"""

  dx = np.float(Lx / (Nx))
  dy = np.float(Ly / (Ny))
  dt = np.float(dx / (2 * c))

  """ defining variable for convenience """

  dt_by_dx = dt / (dx)
  dt_by_dy = dt / (dy)

  """  Defining index grid for updating the fields  """

  I, J = np.meshgrid(range(ghost_cells, x_number_of_points-ghost_cells), \
                      range(ghost_cells, y_number_of_points-ghost_cells)\
                    )

  """  Updating the Electric field  """

  Bz_in_function[I, J] = Bz_in_function[I, J] - ((dt_by_dx * (Ey_in_function[I, J] - Ey_in_function[I, J - 1]))\
                                                  - (dt_by_dy * (Ex_in_function[I, J] - Ex_in_function[I - 1, J]))\
                                                )

  # dBz/dt = - ( dEy/dx - dEx/dy )

  """  Implementing periodic boundary conditions using ghost cells  """

  Bz_in_function = periodic(Bz_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  """  Updating the Magnetic fields   """

  Ex_in_function[I, J] = Ex_in_function[I, J] + (dt_by_dy * (Bz_in_function[I + 1, J] - Bz_in_function[I, J]))

  # dEx/dt = + dBz/dy

  Ey_in_function[I, J] = Ey_in_function[I, J] - (dt_by_dx * (Bz_in_function[I, J + 1] - Bz_in_function[I, J]))

  # dEy/dt = - dBz/dx

  """  Implementing periodic boundary conditions using ghost cells  """

  Ex_in_function = periodic(Ex_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  Ey_in_function = periodic(Ey_in_function, x_number_of_points, y_number_of_points, ghost_cells)

  return Bz_in_function, Ex_in_function, Ey_in_function

"""-------------------------------------------------End--of--Mode--2-------------------------------------------------"""