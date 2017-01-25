import numpy as np
import pylab as pl
from FDTD import fdtd

pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'
pl.rcParams['text.latex.unicode'] = True
pl.rcParams['xtick.major.size'] = 8     
pl.rcParams['xtick.minor.size'] = 4     
pl.rcParams['xtick.major.pad']  = 8     
pl.rcParams['xtick.minor.pad']  = 8     
pl.rcParams['xtick.color']      = 'k'     
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'    
pl.rcParams['image.aspect']  = 'equal'    
pl.rcParams['ytick.major.size'] = 8     
pl.rcParams['ytick.minor.size'] = 4     
pl.rcParams['ytick.major.pad']  = 8     
pl.rcParams['ytick.minor.pad']  = 8     
pl.rcParams['ytick.color']      = 'k'     
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'    

""" User defined function for convinience to find sum of absolute values of all the elements in a 2D matrix"""

def sumsum(a):
    return sum(sum(abs(a) ) )

spread = 0.1
ghost_cells = 1


def error(a, b):

  """ Number of divisions in the physical domain"""
  
  Nx = a
  Ny = b

  """ Speed of light"""

  c = 1
  
  """ Size of domain """
  Lx = 1
  Ly = 1


  """ Setting division size and time steps"""

  dx = np.float(Lx / (Nx))
  dy = np.float(Ly / (Ny))
  dt = np.float(dx / (2 * c))

  """ Setting max number of iterations"""

  max_iterations = np.int(2 / (dt))

  """ defining variable for convinience """

  dt_by_dx = dt / (dx)
  dt_by_dy = dt / (dy)

  """ Setting the spatial physical grids """

  x_center = np.linspace(-dx, 1 + dx, Nx + 1 + 2*ghost_cells, endpoint=True)
  y_center = np.linspace(-dy, 1 + dy, Ny + 1 + 2*ghost_cells, endpoint=True)

  """ Setting the offset spatial grids """

  x_right = np.linspace(-ghost_cells*dx/2, 1 + (2*ghost_cells+1)*dx/2, Nx + 1 + 2*ghost_cells, endpoint=True)
  y_top   = np.linspace(-ghost_cells*dy/2, 1 + (2*ghost_cells+1)*dy/2, Ny + 1 + 2*ghost_cells, endpoint=True)

  """ Writing the spatial grids as a two dimension matrix for vectorization purposes """

  X_center_physical, Y_center_physical = np.meshgrid(x_center[ghost_cells:-ghost_cells], \
                                                  y_center[ghost_cells:-ghost_cells]\
                                                )

  """ Writing the offset spatial grids as a two dimension matrix for vectorization purposes """

  X_right_physical, Y_top_physical     = np.meshgrid(x_right[ghost_cells:-ghost_cells], \
                                                  y_top[ghost_cells:-ghost_cells]\
                                                  )


  I, J = np.meshgrid(range(ghost_cells, len(x_center)-ghost_cells), \
                        range(ghost_cells, len(y_center)-ghost_cells)\
                      )

  """ Initializing the Fields """

  Ez = np.zeros((len(x_center), len(y_center)), dtype=np.float)
  Bx = np.zeros((len(x_center), len(  y_top )), dtype=np.float)
  By = np.zeros((len(x_right ), len(y_center)), dtype=np.float)

  Bz = np.zeros((len(x_center), len(y_center)), dtype=np.float)
  Ex = np.zeros((len(x_center), len(  y_top )), dtype=np.float)
  Ey = np.zeros((len(x_right ), len(y_center)), dtype=np.float)

  Ez_error = 0
  Bx_error = 0
  By_error = 0
    
  Bz_error = 0
  Ex_error = 0
  Ey_error = 0

  """  Setting initial conditions for the fields  """

  Ez[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = np.exp(-((X_center_physical-0.5)**2)/(2*spread**2))

  Bx[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = 0

  By[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = np.exp(-((X_right_physical-0.5)**2)/(2*spread**2))



  Bz[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = 2*np.exp(-((X_right_physical-0.5)**2)/(spread**2))

  Ex[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = 0

  Ey[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] = 2*np.exp(-((X_center_physical-0.5)**2)/(spread**2))
        


  """  Enforcing periodic boundary conditions using ghost cells for all the fields """


  Ez[0, :]                 = Ez[len(y_center) - 2 - ghost_cells, :].copy()
  Ez[:, 0]                 = Ez[:, len(x_center) - 2 - ghost_cells].copy()
  Ez[len(y_center) - 1, :] = Ez[ghost_cells + 1, :].copy()
  Ez[:, len(x_center) - 1] = Ez[:, ghost_cells + 1].copy()

  Bx[0, :]                 = Bx[len(y_top) - 2 - ghost_cells, :].copy()
  Bx[:, 0]                 = Bx[:, len(x_center) - 2 - ghost_cells].copy()
  Bx[len(y_top) - 1, :]    = Bx[ghost_cells + 1, :].copy()
  Bx[:, len(x_center) - 1] = Bx[:, ghost_cells + 1].copy()

  By[0, :]                 = By[len(y_center) - 2 - ghost_cells, :].copy()
  By[:, 0]                 = By[:, len(x_right) - 2 - ghost_cells].copy()
  By[len(y_center) - 1, :] = By[ghost_cells + 1, :].copy()
  By[:, len(x_right) - 1]  = By[:, ghost_cells + 1].copy()
  

  Bz[0, :]                 = Bz[len(y_center) - 2 - ghost_cells, :].copy()
  Bz[:, 0]                 = Bz[:, len(x_center) - 2 - ghost_cells].copy()
  Bz[len(y_center) - 1, :] = Bz[ghost_cells + 1, :].copy()
  Bz[:, len(x_center) - 1] = Bz[:, ghost_cells + 1].copy()

  Ex[0, :]                 = Ex[len(y_top) - 2 - ghost_cells, :].copy()
  Ex[:, 0]                 = Ex[:, len(x_center) - 2 - ghost_cells].copy()
  Ex[len(y_top) - 1, :]    = Ex[ghost_cells + 1, :].copy()
  Ex[:, len(x_center) - 1] = Ex[:, ghost_cells + 1].copy()

  Ey[0, :]                 = Ey[len(y_center) - 2 - ghost_cells, :].copy()
  Ey[:, 0]                 = Ey[:, len(x_right) - 2 - ghost_cells].copy()
  Ey[len(y_center) - 1, :] = Ey[ghost_cells + 1, :].copy()
  Ey[:, len(x_right) - 1]  = Ey[:, ghost_cells + 1].copy()

  Jx = 0
  Jy = 0
  Jz = 0

  """  Starting the solver """

  for time_index in range(max_iterations):
    print('N = ',a,'time = ', time_index)

    Ex, Ey, Ez, Bx, By, Bz = fdtd(Ex, Ey, Ez, Bx, By, Bz, c, Lx, Ly, ghost_cells, Jx, Jy, Jz)


    if(time_index==max_iterations - 1):
      
      Ez_error = sumsum(Ez[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                         np.exp(-((X_center_physical-0.5)**2)/(2*spread**2))\
                       )/(a*b)
      
      Bx_error = sumsum(Bx[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - 0)/(a*b)
      By_error = sumsum(By[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                        np.exp(-((X_right_physical-0.5)**2)/(2*spread**2))\
                       )/(a*b)

      Bz_error = sumsum(Bz[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                        2*np.exp(-((X_right_physical-0.5)**2)/(spread**2))\
                       )/(a*b)
      Ex_error = sumsum(Ex[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - 0)/(a*b)
      
      Ey_error = sumsum(Ey[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] - \
                        2*np.exp(-((X_center_physical-0.5)**2)/(spread**2))\
                       )/(a*b)
      
      return Ez_error, Bx_error, By_error, Bz_error, Ex_error, Ey_error

    
error = np.vectorize(error)


N = np.array( [32, 64, 128, 256] )
ErrorNEz = np.zeros(len(N),dtype = np.float)
ErrorNBx = np.zeros(len(N),dtype = np.float)
ErrorNBy = np.zeros(len(N),dtype = np.float)

ErrorNBz = np.zeros(len(N),dtype = np.float)
ErrorNEx = np.zeros(len(N),dtype = np.float)
ErrorNEy = np.zeros(len(N),dtype = np.float)

ErrorNEz, ErrorNBx, ErrorNBy, ErrorNBz,  ErrorNEx, ErrorNEy,= error(N, N)

pl.loglog(N,ErrorNEz,'-o',label = '$E_z$ ' )
pl.legend()
pl.loglog(N,ErrorNBx,'-o', label = '$B_x$ ' )
pl.legend()
pl.loglog(N,ErrorNBy,'-o', label = '$B_y$ ' )
pl.legend()
pl.loglog(N,ErrorNBz, label = '$B_z$ ' )
pl.legend()
pl.loglog(N,ErrorNEx, label = '$E_x$ ' )
pl.legend()
pl.loglog(N,ErrorNEy, label = '$E_y$ ' )
pl.legend()
pl.loglog(N,1.5*(N**-1.999),'--',color = 'black',lw = 2,label = ' $O(N^{-2})$ ')
pl.legend()
pl.title('$\mathrm{Convergence\; plot}$ ')
pl.xlabel('$\mathrm{N}$')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
pl.show()
pl.clf()        
