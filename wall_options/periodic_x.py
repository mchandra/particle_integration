from simulation_parameters import *
from modules import *

def wall_x(sol):

  x_coordinates = sol[0:no_of_particles]
  wall_x_left   = np.where(x_coordinates<left_boundary)
  wall_x_right  = np.where(x_coordinates>right_boundary)
  
  x_coordinates[wall_x_left[0]]  = x_coordinates[wall_x_left[0]]  + length_box_x
  x_coordinates[wall_x_right[0]] = x_coordinates[wall_x_right[0]] - length_box_x
  
  sol[wall_x_left[0]]  = x_coordinates[wall_x_left[0]]
  sol[wall_x_right[0]] = x_coordinates[wall_x_right[0]]
  return(sol)