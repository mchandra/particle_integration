from ../simulation_parameters import *
from ../modules import *

def wall_y(sol):

  y_coordinates  = sol[no_of_particles:2*no_of_particles]

  if(simulation_dimension == 2): 
    velocity_x = sol[2*no_of_particles:3*no_of_particles]
    velocity_y = sol[3*no_of_particles:4*no_of_particles]

  else:
    velocity_x = sol[3*no_of_particles:4*no_of_particles]
    velocity_y = sol[4*no_of_particles:5*no_of_particles]
    velocity_z = sol[5*no_of_particles:6*no_of_particles]

  collided_top = np.where(y_coordinates>top_boundary)
  collided_bot = np.where(y_coordinates<bottom_boundary)

  # Random variables used to assign the new velocities to the particles:
  R1 = np.random.rand(collided_top[0].size)
  R2 = np.random.rand(collided_top[0].size)

  y_coordinates[collided_top[0]] = top_boundary

  velocity_x[collided_top[0]] = np.sqrt(2*T_top_wall)*erfinv(2*R1-1)    
  velocity_y[collided_top[0]] = np.sqrt(-2*T_top_wall*np.log(R2))*(-1)

  if(simulation_dimension == 3):
    R3                          = np.random.rand(collided_top[0].size)
    velocity_z[collided_top[0]] = np.sqrt(2*T_top_wall)*erfinv(2*R3-1)

  if(simulation_dimension == 2):
    sol[collided_top[0]+2*no_of_particles] = velocity_x[collided_top[0]]
    sol[collided_top[0]+3*no_of_particles] = velocity_y[collided_top[0]]
    
  else:
    sol[collided_top[0]+3*no_of_particles] = velocity_x[collided_top[0]]
    sol[collided_top[0]+4*no_of_particles] = velocity_y[collided_top[0]]
    sol[collided_top[0]+5*no_of_particles] = velocity_z[collided_top[0]]
    
  sol[collided_top[0]+no_of_particles] = y_coordinates[collided_top[0]]

  R1 = np.random.rand(collided_bot[0].size)
  R2 = np.random.rand(collided_bot[0].size)

  y_coordinates[collided_bot[0]] = bottom_boundary

  velocity_x[collided_bot[0]] = np.sqrt(2*T_bot_wall)*erfinv(2*R1-1)    
  velocity_y[collided_bot[0]] = np.sqrt(-2*T_bot_wall*np.log(R2))

  if(simulation_dimension == 3):
    R3                          = np.random.rand(collided_bot[0].size)
    velocity_z[collided_bot[0]] = np.sqrt(2*T_bot_wall)*erfinv(2*R3-1)

  if(simulation_dimension == 2):
    sol[collided_bot[0]+2*no_of_particles] = velocity_x[collided_bot[0]]
    sol[collided_bot[0]+3*no_of_particles] = velocity_y[collided_bot[0]]
    
  else:
    sol[collided_bot[0]+3*no_of_particles] = velocity_x[collided_bot[0]]
    sol[collided_bot[0]+4*no_of_particles] = velocity_y[collided_bot[0]]
    sol[collided_bot[0]+5*no_of_particles] = velocity_z[collided_bot[0]]
    
  sol[collided_bot[0]+no_of_particles] = y_coordinates[collided_bot[0]]
  return(sol)