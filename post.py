from params import *

""" Reading data from generated files"""

# Set the time_range to the number of solution files generated
# varies from case to case
time_range = 400


# Change Nx and Ny to values chosen in field diagnostics.py

Nx = 100
Ny = 100

# Setting dx and dy from Nx and Ny

dx = Lx/Nx
dy = Ly/Ny

x_center = np.linspace(-ghost_cells*dx, Lx + ghost_cells*dx, Nx + 1 + 2 * ghost_cells, endpoint=True)
y_center = np.linspace(-ghost_cells*dy, Ly + ghost_cells*dy, Ny + 1 + 2 * ghost_cells, endpoint=True)

""" Setting the offset spatial grids """

x_right = np.linspace(-ghost_cells * dx / 2, Lx + (2 * ghost_cells + 1) * dx / 2, Nx + 1 + 2 * ghost_cells,\
                        endpoint=True\
                      )

y_top = np.linspace(-ghost_cells * dy / 2, Ly + (2 * ghost_cells + 1) * dy / 2, Ny + 1 + 2 * ghost_cells,\
                      endpoint=True\
                    )



for time_index in range(time_range):

  """ Scripts down below for reading data. Edit these to read the data required"""

  print('post processing for time_index = ', time_index)
  h5f = h5py.File('Ex/solution_'+str(time_index)+'.h5', 'r')
  Ex = h5f['Ex/solution_dataset_'+str(time_index)][:]
  h5f.close()

  h5f = h5py.File('Ey/solution_'+str(time_index)+'.h5', 'r')
  Ey = h5f['Ey/solution_dataset_'+str(time_index)][:]
  h5f.close()

  h5f = h5py.File('Ez/solution_'+str(time_index)+'.h5', 'r')
  Ez = h5f['Ez/solution_dataset_'+str(time_index)][:]
  h5f.close()

  h5f = h5py.File('Bx/solution_'+str(time_index)+'.h5', 'r')
  Bx = h5f['Bx/solution_dataset_'+str(time_index)][:]
  h5f.close()

  h5f = h5py.File('By/solution_'+str(time_index)+'.h5', 'r')
  By = h5f['By/solution_dataset_'+str(time_index)][:]
  h5f.close()

  h5f = h5py.File('Bz/solution_'+str(time_index)+'.h5', 'r')
  Bz = h5f['Bz/solution_dataset_'+str(time_index)][:]
  h5f.close()

  ##movie plots

  # Scripts down below to make movies out of the generated data.
  # Make sure indentation is aligned with the above set of lines
  # By default Ez plots are being saved.

  # script for divergence down below

  #pl.figure(figsize = (13,10) )
  #pl.contourf(x_center[ghost_cells:-ghost_cells], y_center[ghost_cells:-ghost_cells] ,\
              #div_B[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] ,100\
            #)
  #pl.title( r'$\nabla \cdot \mathbf{B}$')
  #pl.xlabel('$x$ ')
  #pl.ylabel('$y$')
  #pl.colorbar()
  #pl.savefig('div/point_mass' + '%04d'%time_index + '.png')  # save it in the folder desired
  ## pl.show()
  #pl.clf()
  #pl.close()



  pl.figure(figsize = (13,10) )
  pl.contourf(x_center[ghost_cells:-ghost_cells], y_center[ghost_cells:-ghost_cells],\
              Ez[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] ,100\
             )
  pl.title(' $E_z(x, y)$')
  pl.xlabel('$x$ ')
  pl.ylabel('$y$')
  pl.colorbar()
  pl.savefig('Ez/point_mass' + '%04d'%time_index + '.png')
  pl.clf()
  pl.close()

  # script for plotting Bz from mode 2 down below:

  #pl.figure(figsize = (13,10) )
  #pl.contourf(x_right[ghost_cells:-ghost_cells], y_top[ghost_cells:-ghost_cells] ,\
              #Bz[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells] ,100\
            #)
  #pl.title(' $B_z(x, y)$')
  #pl.xlabel('$x$ ')
  #pl.ylabel('$y$')
  #pl.colorbar()
  #pl.savefig('Bz/point_mass' + '%04d'%time_index + '.png')
  #pl.clf()
  #pl.close()
