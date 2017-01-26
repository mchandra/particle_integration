from params import *

""" Error function """


def error(a):
  print('Computing error for N = ', a)
  """ Getting the two dimension matrix for initializing the fields """

  Nx = a  # number of zones not points
  Ny = a  # number of zones not points

  """ Length of each zone along x and y """

  dx = np.float(Lx / (Nx))
  dy = np.float(Ly / (Ny))

  x_center = np.linspace(-ghost_cells*dx, 1 + ghost_cells*dx, Nx + 1 + 2 * ghost_cells, endpoint=True)
  y_center = np.linspace(-ghost_cells*dy, 1 + ghost_cells*dy, Ny + 1 + 2 * ghost_cells, endpoint=True)

  """ Initializing the field variables """

  Ez = np.zeros(((len(x_center)), (len(y_center))), dtype=np.float)
  Bx = np.zeros(((len(x_center)), (len(y_center))), dtype=np.float)
  By = np.zeros(((len(x_center)), (len(y_center))), dtype=np.float)


  x_right = np.linspace(-ghost_cells * dx / 2, 1 + (2 * ghost_cells + 1) * dx / 2, Nx + 1 + 2 * ghost_cells,\
                          endpoint=True\
                        )

  y_top = np.linspace(-ghost_cells * dy / 2, 1 + (2 * ghost_cells + 1) * dy / 2, Ny + 1 + 2 * ghost_cells,\
                        endpoint=True\
                      )

  """ Getting the two dimension matrix for initializing the fields """

  X_center_physical, Y_center_physical = np.meshgrid(x_center[ghost_cells:-ghost_cells], \
                                                      y_center[ghost_cells:-ghost_cells] \
                                                    )

  X_right_physical, Y_top_physical = np.meshgrid(x_right[ghost_cells:-ghost_cells], \
                                                  y_top[ghost_cells:-ghost_cells] \
                                                )

  """ [-ghostcells:ghostcells] selects the points located in the physical domain excluding the ghost cells """

  """ Assigning Field values to the physical physical domain """
  Ez[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = initial_fields(X_center_physical, \
                                                                          Y_center_physical
                                                                         )

  Bx[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = initial_fields(X_center_physical, \
                                                                          Y_top_physical \
                                                                         )

  By[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = initial_fields(X_right_physical, \
                                                                          Y_center_physical \
                                                                         )

  """ Implementing Periodic Boundary conditions using ghost cells """

  Ez = periodic(Ez, len(x_center), len(y_center), ghost_cells )
  Bx = periodic(Bx, len(x_center), len(y_top), ghost_cells    )
  By = periodic(By, len(x_right), len(y_center), ghost_cells  )

  """ Selecting a number of test points for testing error """

  number_random_points = 100

  x_random = np.random.rand(number_random_points)
  y_random = np.random.rand(number_random_points)

  """ Selecting a number of test points for testing error """

  Ez_at_random = np.zeros(number_random_points)
  Bx_at_random = np.zeros(number_random_points)
  By_at_random = np.zeros(number_random_points)

  """ Calculating Interpolated values at the randomly selected points """

  Ez_at_random = np.array(bilinear_interpolate(x=[x_random], y=[y_random], x_grid=x_center, \
                                                y_grid=y_center, F=Ez), ghost_cells \
                         )

  Bx_at_random = np.array(bilinear_interpolate(x=[x_random], y=[y_random], x_grid=x_center, \
                                                y_grid=y_top, F=Bx), ghost_cells \
                          )

  By_at_random = np.array(bilinear_interpolate(x=[x_random], y=[y_random], x_grid=x_right, \
                                                y_grid=y_center, F=By), ghost_cells \
                          )

  Ez_error = sum(sum(abs(Ez_at_random - [initial_fields(x_random, y_random)]))) / number_random_points
  Bx_error = sum(sum(abs(Bx_at_random - [initial_fields(x_random, y_random)]))) / number_random_points
  By_error = sum(sum(abs(By_at_random - [initial_fields(x_random, y_random)]))) / number_random_points

  return np.array(Ez_error), np.array(Bx_error), np.array(By_error)


""" Vectorizing the output of the error function"""

error = np.vectorize(error)

""" Choosing test grid densities """

# N = np.array([32, 64, 128, 256, 512, 1024])

N = np.arange(100, 2000, 50)

""" Computing error at the corresponding grid densities """

error_N_Ez, error_N_Bx, error_N_By = error(N)

""" Plotting error vs grid density """

pl.loglog(N, error_N_Ez, '-o', lw=3, label='$E_z$ ')
pl.legend()
pl.loglog(N, error_N_Bx, '-o', lw=3, label='$B_x$ ')
pl.legend()
pl.loglog(N, error_N_By, '-o', lw=5, label='$B_y$ ')
pl.legend()
pl.loglog(N, 150 * (N ** -1.999), '--', color='black', lw=2, label=' $O(N^{-2})$ ')
pl.legend()
pl.title('$\mathrm{Convergence\; plot}$ ')
pl.xlabel('$\mathrm{N}$')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
pl.show()
pl.clf()
