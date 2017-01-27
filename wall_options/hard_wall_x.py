from ../simulation_parameters import *
from ../modules import *

def wall_x(sol):

  x_coordinates  = sol[0:no_of_particles]

  if(simulation_dimension == 2): 
    velocity_x = sol[2*no_of_particles:3*no_of_particles]
  
  else:
    velocity_x = sol[3*no_of_particles:4*no_of_particles]
  
  collided_right = np.where(x_coordinates>right_boundary)
  collided_left  = np.where(x_coordinates<left_boundary)

  x_coordinates[collided_right[0]] = right_boundary
  velocity_x[collided_right[0]]    = velocity_x[collided_right[0]]*(-1)    

  if(simulation_dimension == 2):
    sol[collided_right[0]+2*no_of_particles] = velocity_x[collided_right[0]]
    
  else:
    sol[collided_right[0]+3*no_of_particles] = velocity_x[collided_right[0]]
    
  sol[collided_right[0]] = x_coordinates[collided_right[0]]
  
  x_coordinates[collided_left[0]] = left_boundary
  velocity_x[collided_left[0]]    = velocity_x[collided_left[0]]*(-1)    

  if(simulation_dimension == 2):
    sol[collided_left[0]+2*no_of_particles] = velocity_x[collided_left[0]]
    
  else:
    sol[collided_left[0]+3*no_of_particles] = velocity_x[collided_left[0]]
    
  sol[collided_left[0]] = x_coordinates[collided_left[0]]
  return(sol)