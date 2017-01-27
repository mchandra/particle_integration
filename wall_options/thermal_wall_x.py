from simulation_parameters import *
from modules import *

def wall_x(sol):

  x_coordinates  = sol[0:no_of_particles]

  if(simulation_dimension == 2): 
    velocity_x = sol[2*no_of_particles:3*no_of_particles]
    velocity_y = sol[3*no_of_particles:4*no_of_particles]

  else:
    velocity_x = sol[3*no_of_particles:4*no_of_particles]
    velocity_y = sol[4*no_of_particles:5*no_of_particles]
    velocity_z = sol[5*no_of_particles:6*no_of_particles]

  collided_right = np.where(x_coordinates>right_boundary)
  collided_left  = np.where(x_coordinates<left_boundary)

  # Random variables used to assign the new velocities to the particles:
  R1 = np.random.rand(collided_right[0].size)
  R2 = np.random.rand(collided_right[0].size)

  x_coordinates[collided_right[0]] = right_boundary

  velocity_x[collided_right[0]] = np.sqrt(-2*T_right_wall*np.log(R1))*(-1)    
  velocity_y[collided_right[0]] = np.sqrt(2*T_right_wall)*erfinv(2*R2-1)

  if(simulation_dimension == 3):
    R3                            = np.random.rand(collided_right[0].size)
    velocity_z[collided_right[0]] = np.sqrt(2*T_right_wall)*erfinv(2*R3-1)

  if(simulation_dimension == 2):
    sol[collided_right[0]+2*no_of_particles] = velocity_x[collided_right[0]]
    sol[collided_right[0]+3*no_of_particles] = velocity_y[collided_right[0]]
    
  else:
    sol[collided_right[0]+3*no_of_particles] = velocity_x[collided_right[0]]
    sol[collided_right[0]+4*no_of_particles] = velocity_y[collided_right[0]]
    sol[collided_right[0]+5*no_of_particles] = velocity_z[collided_right[0]]
    
  sol[collided_right[0]] = x_coordinates[collided_right[0]]

  R1 = np.random.rand(collided_left[0].size)
  R2 = np.random.rand(collided_left[0].size)

  x_coordinates[collided_left[0]] = left_boundary

  velocity_x[collided_left[0]] = np.sqrt(-2*T_left_wall*np.log(R1))   
  velocity_y[collided_left[0]] = np.sqrt(2*T_left_wall)*erfinv(2*R2-1)

  if(simulation_dimension == 3):
    R3                           = np.random.rand(collided_left[0].size)
    velocity_z[collided_left[0]] = np.sqrt(2*T_left_wall)*erfinv(2*R3-1)

  if(simulation_dimension == 2):
    sol[collided_left[0]+2*no_of_particles] = velocity_x[collided_left[0]]
    sol[collided_left[0]+3*no_of_particles] = velocity_y[collided_left[0]]
    
  else:
    sol[collided_left[0]+3*no_of_particles] = velocity_x[collided_left[0]]
    sol[collided_left[0]+4*no_of_particles] = velocity_y[collided_left[0]]
    sol[collided_left[0]+5*no_of_particles] = velocity_z[collided_left[0]]
    
  sol[collided_left[0]] = x_coordinates[collided_left[0]]
  return(sol)