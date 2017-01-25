import numpy as np
import pylab as pl
import numpy.linalg as la
from Interpolator import bilinear_interpolate

pl.rcParams['figure.figsize'] = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family'] = 'serif'
pl.rcParams['font.weight'] = 'bold'
pl.rcParams['font.size'] = 20
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex'] = True
pl.rcParams['axes.linewidth'] = 1.5
pl.rcParams['axes.titlesize'] = 'medium'
pl.rcParams['axes.labelsize'] = 'medium'
pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.major.pad'] = 8
pl.rcParams['xtick.minor.pad'] = 8
pl.rcParams['xtick.color'] = 'k'
pl.rcParams['xtick.labelsize'] = 'medium'
pl.rcParams['xtick.direction'] = 'in'
pl.rcParams['ytick.major.size'] = 8
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.major.pad'] = 8
pl.rcParams['ytick.minor.pad'] = 8
pl.rcParams['ytick.color'] = 'k'
pl.rcParams['ytick.labelsize'] = 'medium'
pl.rcParams['ytick.direction'] = 'in'

""" Number of Ghost cells """

ghost_cells = 1

""" Function for setting initial conditions for the fields """

def initial_fields(x, y):
  function_value = np.sin(2 * np.pi * x * y) * np.cos(2 * np.pi * x * y)

  return function_value


""" Error function """


def error(a):
  print('Computing error for N = ', a)
  """ Getting the two dimension matrix for initializing the fields """

  nx = a  # number of zones not points
  ny = a  # number of zones not points

  """ Length of each zone along x and y """

  dx = np.float(1 / (nx))
  dy = np.float(1 / (ny))

  x_center = np.linspace(-dx, 1 + dx, nx + 3, endpoint=True)
  y_center = np.linspace(-dy, 1 + dy, ny + 3, endpoint=True)

  """ Initializing the field variables """

  Ez = np.zeros(((len(x_center)), (len(y_center))), dtype=np.float)
  Bx = np.zeros(((len(x_center)), (len(y_center))), dtype=np.float)
  By = np.zeros(((len(x_center)), (len(y_center))), dtype=np.float)

  x_right = np.linspace(-dx / 2, 1 + 3 * dx / 2, nx + 3, endpoint=True)
  y_top = np.linspace(-dy / 2, 1 + 3 * dy / 2, ny + 3, endpoint=True)

  """ Getting the two dimension matrix for initializing the fields """

  X_center_physical, Y_center_physical = np.meshgrid(x_center[ghost_cells:-ghost_cells],\
                                                      y_center[ghost_cells:-ghost_cells]\
                                                    )
  
  X_right_physical, Y_top_physical = np.meshgrid(x_right[ghost_cells:-ghost_cells],\
                                                 y_top[ghost_cells:-ghost_cells]\
                                                )

  """ [-ghostcells:ghostcells] selects the points located in the physical domain excluding the ghost cells """

  """ Assigning Field values to the physical physical domain """
  Ez[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = initial_fields(X_center_physical,\
                                                                          Y_center_physical
                                                                         )
  
  Bx[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = initial_fields(X_center_physical,\
                                                                          Y_top_physical\
                                                                         )
  
  By[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = initial_fields(X_right_physical,\
                                                                          Y_center_physical\
                                                                         )

  """ Implementing Periodic Boundary conditions using ghost cells """

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

  """ Selecting a number of test points for testing error """

  number_random_points = 100

  x_random = np.random.rand(number_random_points)
  y_random = np.random.rand(number_random_points)

  """ Selecting a number of test points for testing error """

  Ez_at_random = np.zeros(number_random_points)
  Bx_at_random = np.zeros(number_random_points)
  By_at_random = np.zeros(number_random_points)

  """ Calculating Interpolated values at the randomly selected points """

  Ez_at_random = np.array(bilinear_interpolate(x=[x_random], y=[y_random], x_grid=x_center,\
                                               y_grid=y_center, F=Ez)\
                                              )
  
  Bx_at_random = np.array(bilinear_interpolate(x=[x_random], y=[y_random], x_grid=x_center,\
                                               y_grid=y_top, F=Bx)\
                                              )
  
  By_at_random = np.array(bilinear_interpolate(x=[x_random], y=[y_random], x_grid=x_right,\
                                              y_grid=y_center, F=By)\
                                              )

  Ez_error = sum(sum(abs(Ez_at_random - [initial_fields(x_random, y_random)]))) / number_random_points
  Bx_error = sum(sum(abs(Bx_at_random - [initial_fields(x_random, y_random)]))) / number_random_points
  By_error = sum(sum(abs(By_at_random - [initial_fields(x_random, y_random)]))) / number_random_points

  return np.array(Ez_error), np.array(Bx_error), np.array(By_error)


""" Vectorizing the output of the error function"""

error = np.vectorize(error)

""" Choosing test grid densities """

# N = np.array([32, 64, 128, 256, 512, 1024])

N = np.arange(100, 4000, 100)

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
