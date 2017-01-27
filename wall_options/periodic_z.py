from simulation_parameters import *
from modules import *

def wall_z(sol):

  z_coordinates = sol[2*no_of_particles:3*no_of_particles]
  wall_z_back   = np.where(z_coordinates<back_boundary)
  wall_z_front  = np.where(z_coordinates>front_boundary)
  
  z_coordinates[wall_z_back[0]]  = z_coordinates[wall_z_back[0]]  + length_box_z
  z_coordinates[wall_z_front[0]] = z_coordinates[wall_z_front[0]] - length_box_z
  
  sol[2*no_of_particles+wall_z_back[0]]  = z_coordinates[wall_z_back[0]]
  sol[2*no_of_particles+wall_z_front[0]] = z_coordinates[wall_z_front[0]]
  return(sol)