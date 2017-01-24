import numpy as np
import h5py
import pylab as pl

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

"""  Setting number of ghost cells  """

spread = 0.1
ghost_cells = 1


def two_d_fdtd(mode, Nx, Ny, t):

    if mode==1:
        return mode1_fdtd(Nx, Ny, t)
    else:
        return mode2_fdtd(Nx, Ny, t)





""" Equations for mode 1 fdtd (variation along x and y)"""

# dEz/dt = dBy/dx - dBx/dy
# dBx/dt = -dEz/dy
# dBy/dt = +dEz/dx
# div_B = dBx/dx + dBy/dy


"""Particle grid coinciding with  the spatial grid for the Electric field """

def mode1_fdtd(a, b, time_in_seconds):

    """ Number of divisions in the physical domain"""
    Nx = a
    Ny = b

    """ Speed of light"""

    c = 1

    """ Setting division size and time steps"""

    dx = np.float(1 / (Nx))
    dy = np.float(1 / (Ny))
    dt = np.float(dx / (2 * c))

    """ Setting max number of iterations"""

    max_iterations = np.int(time_in_seconds / (dt))

    """ defining variable for convinience """

    dt_by_dx = dt / (dx)
    dt_by_dy = dt / (dy)

    """ Setting the spatial physical grids """

    x_center = np.linspace(-dx, 1 + dx, Nx + 3, endpoint=True)
    y_center = np.linspace(-dy, 1 + dy, Ny + 3, endpoint=True)

    x_right = np.linspace(-ghost_cells*dx/2, 1 + (2*ghost_cells+1)*dx/2, Nx + 1 + 2*ghost_cells, endpoint=True)
    y_top   = np.linspace(-ghost_cells*dy/2, 1 + (2*ghost_cells+1)*dy/2, Ny + 1 + 2*ghost_cells, endpoint=True)

    """ Writing the spatial grids as a two dimension matrix for vectorization purposes """

    X_center_physical, Y_center_physical = np.meshgrid(x_center[ghost_cells:-ghost_cells], \
                                                   y_center[ghost_cells:-ghost_cells]
                                                  )

    """ Writing the offset spatial grids as a two dimension matrix for vectorization purposes """

    X_right_physical, Y_top_physical     = np.meshgrid(x_right[ghost_cells:-ghost_cells], \
                                                   y_top[ghost_cells:-ghost_cells]
                                                   )

    """ Initializing the Fields """

    Ez = np.zeros((len(x_center), len(y_center)), dtype=np.float)
    Bx = np.zeros((len(x_center), len(  y_top )), dtype=np.float)
    By = np.zeros((len(x_right ), len(y_center)), dtype=np.float)

    """ Initializing the  divergence matrix (Optional)"""

    div_B = np.zeros((len(x_center), len(  y_top )), dtype=np.float)

    """  Starting the solver """

    for time_index in range(max_iterations):
        print('Mode 1 with time = ', time_index)

        if (time_index == 0):

            """  Setting initial conditions for the fields  """

            Ez[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = np.exp(-( (Y_center_physical - 0.5) ** 2 + \
                                                                              (X_center_physical - 0.5) ** 2) /\
                                                                               (2 * spread ** 2)
                                                                           )

            Bx[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = np.exp(-((Y_top_physical - 0.5) ** 2) /\
                                                                           (2 * spread ** 2)
                                                                           )

            By[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = np.exp(-((X_right_physical - 0.5) ** 2) /\
                                                                           (2 * spread ** 2)
                                                                           )

            """  Implementing periodic boundary conditions using ghost cells for all the fields """

            Ez[0, :]                 = Ez[len(y_center) - 1 - ghost_cells, :].copy()
            Ez[:, 0]                 = Ez[:, len(x_center) - 1 - ghost_cells].copy()
            Ez[len(y_center) - 1, :] = Ez[ghost_cells, :].copy()
            Ez[:, len(x_center) - 1] = Ez[:, ghost_cells].copy()

            Bx[0, :]                 = Bx[len(y_top) - 1 - ghost_cells, :].copy()
            Bx[:, 0]                 = Bx[:, len(x_center) - 1 - ghost_cells].copy()
            Bx[len(y_top) - 1, :]    = Bx[ghost_cells, :].copy()
            Bx[:, len(x_center) - 1] = Bx[:, ghost_cells].copy()

            By[0, :]                 = By[len(y_center) - 1 - ghost_cells, :].copy()
            By[:, 0]                 = By[:, len(x_right) - 1 - ghost_cells].copy()
            By[len(y_center) - 1, :] = By[ghost_cells, :].copy()
            By[:, len(x_right) - 1]  = By[:, ghost_cells].copy()

        """  Updating the Electric field  """

        I, J = np.meshgrid(range(ghost_cells, len(x_center)-ghost_cells), \
                           range(ghost_cells, len(y_center)-ghost_cells)
                           )

        Ez[I, J] = Ez[I, J] + (dt_by_dx * (By[I, J] - By[I, J - 1]))\
                   - (dt_by_dy * (Bx[I, J] - Bx[I - 1, J])
                      )

        """  Implementing periodic boundary conditions using ghost cells  """

        Ez[0, :]                 = Ez[len(y_center) - 1 - ghost_cells, :].copy()
        Ez[:, 0]                 = Ez[:, len(x_center) - 1 - ghost_cells].copy()
        Ez[len(y_center) - 1, :] = Ez[ghost_cells, :].copy()
        Ez[:, len(x_center) - 1] = Ez[:, ghost_cells].copy()

        """  Updating the Magnetic fields   """

        Bx[I, J] = Bx[I, J] - (dt_by_dy * (Ez[I + 1, J] - Ez[I, J]))

        By[I, J] = By[I, J] + (dt_by_dx * (Ez[I, J + 1] - Ez[I, J]))

        """  Implementing periodic boundary conditions using ghost cells  """

        Bx[0, :]                 = Bx[len(y_top) - 1 - ghost_cells, :].copy()
        Bx[:, 0]                 = Bx[:, len(x_center) - 1 - ghost_cells].copy()
        Bx[len(y_top) - 1, :]    = Bx[ghost_cells, :].copy()
        Bx[:, len(x_center) - 1] = Bx[:, ghost_cells].copy()



        By[0, :]                 = By[len(y_center) - 1 - ghost_cells, :].copy()
        By[:, 0]                 = By[:, len(x_right) - 1 - ghost_cells].copy()
        By[len(y_center) - 1, :] = By[ghost_cells, :].copy()
        By[:, len(x_right) - 1]  = By[:, ghost_cells].copy()

        """  Calculating the current divergence (Optional)  """

        div_B[I, J] = (Bx[I, J + 1] - Bx[I, J]) / (dx) + (By[I + 1, J] - By[I, J]) / (dy)

        div_B[0, :]                 = div_B[len(y_center) - 1 - ghost_cells, :].copy()
        div_B[:, 0]                 = div_B[:, len(x_right) - 1 - ghost_cells].copy()
        div_B[len(y_center) - 1, :] = div_B[ghost_cells, :].copy()
        div_B[:, len(x_right) - 1]  = div_B[:, ghost_cells].copy()


"""-------------------------------------------------End--of--Mode--1-------------------------------------------------"""



""" Equations for mode 2 fdtd (variation along x and y)"""

"""Particle grid coinciding with  the spatial grid for the Magnetic field """

# dBz/dt = - ( dEy/dx - dEx/dy )
# dEx/dt = + dBz/dy
# dEy/dt = - dBz/dx
# div_B = dBz/dz


def mode2_fdtd(a, b, time_in_seconds):

    """ Number of divisions in the physical domain"""

    Nx = a
    Ny = b

    """ Speed of light"""

    c = 1

    """ Setting division size and time steps"""

    dx = np.float(1 / (Nx))
    dy = np.float(1 / (Ny))
    dt = np.float(dx / (2 * c))

    """ Setting max number of iterations"""

    max_iterations = np.int(time_in_seconds / (dt))

    """ defining variable for convinience """

    dt_by_dx = dt / (dx)
    dt_by_dy = dt / (dy)

    """ Setting the spatial physical grids """

    x_center = np.linspace(-dx, 1 + dx, Nx + 3, endpoint=True)
    y_center = np.linspace(-dy, 1 + dy, Ny + 3, endpoint=True)

    x_right = np.linspace(-ghost_cells*dx/2, 1 + (2*ghost_cells+1)*dx/2, Nx + 1 + 2*ghost_cells, endpoint=True)
    y_top   = np.linspace(-ghost_cells*dy/2, 1 + (2*ghost_cells+1)*dy/2, Ny + 1 + 2*ghost_cells, endpoint=True)

    """ Writing the spatial grids as a two dimension matrix for vectorization purposes """

    X_center_physical, Y_center_physical = np.meshgrid(x_center[ghost_cells:-ghost_cells], \
                                                   y_center[ghost_cells:-ghost_cells]
                                                  )

    """ Writing the offset spatial grids as a two dimension matrix for vectorization purposes """

    X_right_physical, Y_top_physical     = np.meshgrid(x_right[ghost_cells:-ghost_cells], \
                                                   y_top[ghost_cells:-ghost_cells]
                                                   )

    """ Initializing the Fields """

    Bz = np.zeros((len(x_center), len(y_center)), dtype=np.float)
    Ex = np.zeros((len(x_center), len(  y_top )), dtype=np.float)
    Ey = np.zeros((len(x_right ), len(y_center)), dtype=np.float)

    """  Starting the solver """

    for time_index in range(max_iterations):
        print('Mode 2 with Time = ', time_index)

        if (time_index == 0):

            """  Setting initial conditions for the fields  """

            Bz[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = np.exp(-( (Y_center_physical - 0.5) ** 2 + \
                                                                              (X_center_physical - 0.5) ** 2) /\
                                                                               (2 * spread ** 2)
                                                                           )

            Ex[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = np.exp(-((Y_top_physical - 0.5) ** 2) /\
                                                                           (2 * spread ** 2)
                                                                           )

            Ey[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = np.exp(-((X_right_physical - 0.5) ** 2) /\
                                                                           (2 * spread ** 2)
                                                                           )

            """  Implementing periodic boundary conditions using ghost cells for all the fields """

            Bz[0, :]                 = Bz[len(y_center) - 1 - ghost_cells, :].copy()
            Bz[:, 0]                 = Bz[:, len(x_center) - 1 - ghost_cells].copy()
            Bz[len(y_center) - 1, :] = Bz[ghost_cells, :].copy()
            Bz[:, len(x_center) - 1] = Bz[:, ghost_cells].copy()

            Ex[0, :]                 = Ex[len(y_top) - 1 - ghost_cells, :].copy()
            Ex[:, 0]                 = Ex[:, len(x_center) - 1 - ghost_cells].copy()
            Ex[len(y_top) - 1, :]    = Ex[ghost_cells, :].copy()
            Ex[:, len(x_center) - 1] = Ex[:, ghost_cells].copy()

            Ey[0, :]                 =     Ey[len(y_center) - 1 - ghost_cells, :].copy()
            Ey[:, 0]                 =     Ey[:, len(x_right) - 1 - ghost_cells].copy()
            Ey[len(y_center) - 1, :] =     Ey[ghost_cells, :].copy()
            Ey[:, len(x_right) - 1]  =     Ey[:, ghost_cells].copy()

        """  Updating the Magnetic field  """

        I, J = np.meshgrid(range(ghost_cells, len(x_center)-ghost_cells), \
                           range(ghost_cells, len(y_center)-ghost_cells)
                           )

        Bz[I, J] = Bz[I, J] -((dt_by_dx * (    Ey[I, J] -     Ey[I, J - 1]))\
                               - (dt_by_dy * (Ex[I, J] - Ex[I - 1, J]))\
                               )

        """  Implementing periodic boundary conditions using ghost cells  """

        Bz[0, :]                 = Bz[len(y_center) - 1 - ghost_cells, :].copy()
        Bz[:, 0]                 = Bz[:, len(x_center) - 1 - ghost_cells].copy()
        Bz[len(y_center) - 1, :] = Bz[ghost_cells, :].copy()
        Bz[:, len(x_center) - 1] = Bz[:, ghost_cells].copy()

        """  Updating the Electric fields   """

        Ex[I, J] = Ex[I, J] + (dt_by_dy * (Bz[I + 1, J] - Bz[I, J]))

        Ey[I, J] =     Ey[I, J] - (dt_by_dx * (Bz[I, J + 1] - Bz[I, J]))

        """  Implementing periodic boundary conditions using ghost cells  """

        Ex[0, :]                 = Ex[len(y_top) - 1 - ghost_cells, :].copy()
        Ex[:, 0]                 = Ex[:, len(x_center) - 1 - ghost_cells].copy()
        Ex[len(y_top) - 1, :]    = Ex[ghost_cells, :].copy()
        Ex[:, len(x_center) - 1] = Ex[:, ghost_cells].copy()



        Ey[0, :]                 =     Ey[len(y_center) - 1 - ghost_cells, :].copy()
        Ey[:, 0]                 =     Ey[:, len(x_right) - 1 - ghost_cells].copy()
        Ey[len(y_center) - 1, :] =     Ey[ghost_cells, :].copy()
        Ey[:, len(x_right) - 1]  =     Ey[:, ghost_cells].copy()

"""-------------------------------------------------End--of--Mode--2-------------------------------------------------"""


# Example for running the fdtd code

#f = two_d_fdtd(1,100,100) # for mode 1
#f = two_d_fdtd(2,100,100, 4) # for mode 2
