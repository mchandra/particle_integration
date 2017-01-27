from ../simulation_parameters import *
from ../modules import *

def wall_z(sol):

  z_coordinates = sol[2*no_of_particles:3*no_of_particles]
  velocity_x    = sol[3*no_of_particles:4*no_of_particles]
  velocity_y    = sol[4*no_of_particles:5*no_of_particles]
  velocity_z    = sol[5*no_of_particles:6*no_of_particles]

  collided_front = np.where(z_coordinates>front_boundary)
  collided_back  = np.where(z_coordinates<back_boundary)

  # Random variables used to assign the new velocities to the particles:
  R1 = np.random.rand(collided_front[0].size)
  R2 = np.random.rand(collided_front[0].size)
  R3 = np.random.rand(collided_front[0].size)

  z_coordinates[collided_front[0]] = front_boundary

  velocity_x[collided_front[0]] = np.sqrt(2*T_front_wall)*erfinv(2*R1-1)    
  velocity_y[collided_front[0]] = np.sqrt(2*T_front_wall)*erfinv(2*R2-1)
  velocity_z[collided_front[0]] = np.sqrt(-2*T_front_wall*np.log(R3))*(-1)

  sol[collided_front[0]+3*no_of_particles] = velocity_x[collided_front[0]]
  sol[collided_front[0]+4*no_of_particles] = velocity_y[collided_front[0]]
  sol[collided_front[0]+5*no_of_particles] = velocity_z[collided_front[0]]
    
  sol[collided_front[0]+2*no_of_particles] = z_coordinates[collided_front[0]]

  R1 = np.random.rand(collided_back[0].size)
  R2 = np.random.rand(collided_back[0].size)
  R3 = np.random.rand(collided_back[0].size)

  z_coordinates[collided_back[0]] = back_boundary

  velocity_x[collided_back[0]] = np.sqrt(2*T_back_wall)*erfinv(2*R1-1)    
  velocity_y[collided_back[0]] = np.sqrt(2*T_back_wall)*erfinv(2*R2-1)
  velocity_z[collided_back[0]] = np.sqrt(-2*T_back_wall*np.log(R3))

  sol[collided_back[0]+3*no_of_particles] = velocity_x[collided_back[0]]
  sol[collided_back[0]+4*no_of_particles] = velocity_y[collided_back[0]]
  sol[collided_back[0]+5*no_of_particles] = velocity_z[collided_back[0]]
    
  sol[collided_back[0]+2*no_of_particles] = z_coordinates[collided_back[0]]
  
  return(sol)