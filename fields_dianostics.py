from params import *

def error_convergence(a, b):
  """ Number of divisions in the physical domain"""

  Nx = a
  Ny = b

  """ Setting division size and time steps"""

  dx = np.float(Lx / (Nx))
  dy = np.float(Ly / (Ny))

  dt = np.float(dx / (2 * c))

  """ Setting max number of iterations"""

  max_iterations = np.int(2 / (dt))

  """ Setting the spatial physical grids """

  x_center = np.linspace(-ghost_cells*dx, 1 + ghost_cells*dx, Nx + 1 + 2 * ghost_cells, endpoint=True)
  y_center = np.linspace(-ghost_cells*dy, 1 + ghost_cells*dy, Ny + 1 + 2 * ghost_cells, endpoint=True)

  """ Setting the offset spatial grids """

  x_right = np.linspace(-ghost_cells * dx / 2, 1 + (2 * ghost_cells + 1) * dx / 2, Nx + 1 + 2 * ghost_cells,\
                          endpoint=True\
                        )

  y_top = np.linspace(-ghost_cells * dy / 2, 1 + (2 * ghost_cells + 1) * dy / 2, Ny + 1 + 2 * ghost_cells,\
                        endpoint=True\
                      )

  """ Writing the spatial grids as a two dimension matrix for vectorization purposes """

  X_center_physical, Y_center_physical = np.meshgrid( x_center[ghost_cells:-ghost_cells], \
                                                      y_center[ghost_cells:-ghost_cells] \
                                                    )

  """ Writing the offset spatial grids as a two dimension matrix for vectorization purposes """

  X_right_physical, Y_top_physical = np.meshgrid( x_right[ghost_cells:-ghost_cells], \
                                                  y_top[ghost_cells:-ghost_cells] \
                                                )

  I, J = np.meshgrid( range(ghost_cells, len(x_center) - ghost_cells), \
                      range(ghost_cells, len(y_center) - ghost_cells) \
                    )

  """ Initializing the Fields """

  Ez = np.zeros((len(x_center), len(y_center)), dtype=np.float)
  Bx = np.zeros((len(x_center), len(y_top)   ), dtype=np.float)
  By = np.zeros((len(x_right), len(y_center) ), dtype=np.float)

  Bz = np.zeros((len(x_center), len(y_center)), dtype=np.float)
  Ex = np.zeros((len(x_center), len(y_top)   ), dtype=np.float)
  Ey = np.zeros((len(x_right), len(y_center) ), dtype=np.float)

  Ez_initial = np.zeros((len(x_center), len(y_center)), dtype=np.float)
  Bx_initial = np.zeros((len(x_center), len(y_top)   ), dtype=np.float)
  By_initial = np.zeros((len(x_right), len(y_center) ), dtype=np.float)

  Bz_initial = np.zeros((len(x_center), len(y_center)), dtype=np.float)
  Ex_initial = np.zeros((len(x_center), len(y_top)   ), dtype=np.float)
  Ey_initial = np.zeros((len(x_right), len(y_center) ), dtype=np.float)

  """  Setting initial conditions for the fields in the physical domain """

  # Initialize the fields in the manner desired below:

  Ez[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = gauss1D(X_center_physical)

  Bx[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = 0

  By[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = gauss1D(X_right_physical)

  Bz[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = gauss1D(2*X_right_physical)

  Ex[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = 0

  Ey[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = gauss1D(2*X_center_physical)

  # Saving the initial conditions in different variables

  Ez_initial = Ez.copy()
  Bx_initial = Bx.copy()
  By_initial = By.copy()

  Bz_initial = Bz.copy()
  Ex_initial = Ex.copy()
  Ey_initial = Ey.copy()

  """  Starting the solver """

  for time_index in range(max_iterations):
    print('N = ', a, 'time = ', time_index)

    Ex, Ey, Ez, Bx, By, Bz = fdtd(Ex, Ey, Ez, Bx, By, Bz, c, Lx, Ly, ghost_cells, Jx, Jy, Jz)

    # For arbitrary initial conditions set time_index == timestep # where wave comes back to its initial conditions
    # if it is happening

    if (time_index == max_iterations - 1):
      
      Ez_error = sumsum( Ez[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                         Ez_initial[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells]\
                       ) / (a * b)

      Bx_error = sumsum( Bx[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                         Bx_initial[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells]\
                       ) / (a * b)


      By_error = sumsum( By[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                         By_initial[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells]\
                       ) / (a * b)

      Bz_error = sumsum( Bz[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                         Bz_initial[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells]\
                       ) / (a * b)

      Ex_error = sumsum( Ex[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                         Ex_initial[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells]\
                       ) / (a * b)

      Ey_error =  sumsum( Ey[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                          Ey_initial[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells]\
                       ) / (a * b)

      return Ez_error, Bx_error, By_error, Bz_error, Ex_error, Ey_error


      # Optional Script for making movies ()
      # Add scripts here for resulting in images for variables desired
      # Adjust the scripts according to your needs

      """












      """

error_convergence = np.vectorize(error_convergence)

N = np.array([32, 64, 128, 256])

ErrorNEz, ErrorNBx, ErrorNBy, ErrorNBz, ErrorNEx, ErrorNEy, = error_convergence(N, N)

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
