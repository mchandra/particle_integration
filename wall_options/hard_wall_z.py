from simulation_parameters import *
from modules import *

def wall_z(sol):

  z_coordinates  = sol[2*no_of_particles:3*no_of_particles]
  velocity_z     = sol[5*no_of_particles:6*no_of_particles]

  collided_front = np.where(z_coordinates>front_boundary)
  collided_back  = np.where(z_coordinates<back_boundary)

  z_coordinates[collided_front[0]] = front_boundary
  velocity_z[collided_front[0]]    = velocity_z[collided_front[0]]*(-1)

  sol[collided_front[0]+5*no_of_particles] = velocity_z[collided_front[0]]
  sol[collided_front[0]+2*no_of_particles] = z_coordinates[collided_front[0]]

  z_coordinates[collided_back[0]] = back_boundary
  velocity_z[collided_back[0]]    = velocity_z[collided_back[0]]*(-1)

  sol[collided_back[0]+5*no_of_particles] = velocity_z[collided_back[0]]
  sol[collided_back[0]+2*no_of_particles] = z_coordinates[collided_back[0]]
  
  return(sol)