import numpy as np
import h5py
import pylab as pl

""" Equations for mode 1 fdtd"""

# dEz/dt = dBy/dx - dBx/dy
# dBx/dt = -dEz/dy
# dBy/dt = +dEz/dx
# div_B = dBx/dx + dBy/dy

"""  Setting number of ghost cells  """

spread = 0.1
ghost_cells = 1


def mode1_fdtd(a, b):
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

    max_iterations = np.int(2 / (dt))

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
        print('Time = ', time_index)

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

        """  Calculating the current divergence   """

        div_B[I, J] = (Bx[I, J + 1] - Bx[I, J]) / (dx) + (By[I + 1, J] - By[I, J]) / (dy)

        div_B[0, :]                 = div_B[len(y_center) - 1 - ghost_cells, :].copy()
        div_B[:, 0]                 = div_B[:, len(x_right) - 1 - ghost_cells].copy()
        div_B[len(y_center) - 1, :] = div_B[ghost_cells, :].copy()
        div_B[:, len(x_right) - 1]  = div_B[:, ghost_cells].copy()

