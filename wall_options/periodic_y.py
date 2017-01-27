from simulation_parameters import *
from modules import *

def wall_y(sol):

  y_coordinates = sol[no_of_particles:2*no_of_particles]
  wall_y_bot    = np.where(y_coordinates<bottom_boundary)
  wall_y_top    = np.where(y_coordinates>top_boundary)
  
  y_coordinates[wall_y_bot[0]] = y_coordinates[wall_y_bot[0]] + length_box_y
  y_coordinates[wall_y_top[0]] = y_coordinates[wall_y_top[0]] - length_box_y
  
  sol[no_of_particles+wall_y_bot[0]] = y_coordinates[wall_y_bot[0]]
  sol[no_of_particles+wall_y_top[0]] = y_coordinates[wall_y_top[0]]
  return(sol)