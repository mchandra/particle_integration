from ../simulation_parameters import *
from ../modules import *

def wall_y(sol):

  y_coordinates  = sol[no_of_particles:2*no_of_particles]

  if(simulation_dimension == 2): 
    velocity_y = sol[3*no_of_particles:4*no_of_particles]

  else:
    velocity_y = sol[4*no_of_particles:5*no_of_particles]
  
  collided_top = np.where(y_coordinates>top_boundary)
  collided_bot = np.where(y_coordinates<bottom_boundary)

  y_coordinates[collided_top[0]] = top_boundary
  velocity_y[collided_top[0]]    = velocity_y[collided_top[0]]*(-1)

  if(simulation_dimension == 2):
    sol[collided_top[0]+3*no_of_particles] = velocity_x[collided_top[0]]
    
  else:
    sol[collided_top[0]+4*no_of_particles] = velocity_x[collided_top[0]]
    
  sol[collided_top[0]+no_of_particles] = y_coordinates[collided_top[0]]

  y_coordinates[collided_bot[0]] = bottom_boundary
  velocity_y[collided_bot[0]]    = velocity_y[collided_bot[0]]*(-1)

  if(simulation_dimension == 2):
    sol[collided_bot[0]+3*no_of_particles] = velocity_x[collided_bot[0]]
    
  else:
    sol[collided_bot[0]+4*no_of_particles] = velocity_x[collided_bot[0]]
    
  sol[collided_bot[0]+no_of_particles] = y_coordinates[collided_bot[0]]
  return(sol)