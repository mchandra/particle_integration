import numpy as np
from scipy.special import erfinv
import h5py
import params

""" 
This file contains 3 functions, which define periodic B.C's in 3 directions
Depending upon the choice of the user, periodic boundary conditions may be set
to either of the x,y and z directions.

A periodic B.C means that a particle that encounters such a boundary will result
in the particle reaching the other side of the box. Effectively, this means that
the particle in moving on a ring
"""

"""Here we shall assign values as set in params"""

no_of_particles      = params.no_of_particles
simulation_dimension = params.simulation_dimension
restart_simulation   = params.restart_simulation
choice_integrator    = params.choice_integrator
collision_operator   = params.collision_operator

if(collision_operator == "hardsphere"):
  scattering_distance = params.scattering_distance

elif(collision_operator == "potential-based"):
  potential_steepness     = params.potential_steepness
  potential_amplitude     = params.potential_amplitude
  order_finite_difference = params.order_finite_difference

elif(collision_operator == "montecarlo"):
  x_zones            = params.x_zones
  y_zones            = params.y_zones
  scattered_fraction = params.scattered_fraction

mass_particle      = params.mass_particle
boltzmann_constant = params.boltzmann_constant
T_initial          = params.T_initial
wall_condition_x   = params.wall_condition_x
wall_condition_y   = params.wall_condition_y
wall_condition_z   = params.wall_condition_z

if(wall_condition_x == "thermal"):
  T_left_wall  = params.T_left_wall
  T_right_wall = params.T_right_wall

if(wall_condition_y == "thermal"):
  T_top_wall = params.T_top_wall
  T_bot_wall = params.T_bot_wall

if(wall_condition_z == "thermal"):
  T_front_wall = params.T_front_wall
  T_back_wall  = params.T_back_wall

left_boundary    = params.left_boundary
right_boundary   = params.right_boundary
length_box_x     = params.length_box_x

bottom_boundary  = params.bottom_boundary
top_boundary     = params.top_boundary
length_box_y     = params.length_box_y

back_boundary    = params.back_boundary
front_boundary   = params.front_boundary
length_box_z     = params.length_box_z

#Here we complete import of all the variable from the parameters file

def periodic(coordinates, direction):

  if(direction == 'x'):
    boundary          = right_boundary
    boundary_opposite = left_boundary

  if(direction == 'y'):
    boundary          = top_boundary
    boundary_opposite = bottom_boundary

  if(direction == 'z'):
    boundary          = front_boundary
    boundary_opposite = back_boundary

  collided          = np.where(coordinates > boundary)
  collided_opposite = np.where(coordinates < boundary_opposite)
  
  coordinates[collided[0]]          = coordinates[collided[0]]  - length_box_x
  coordinates[collided_opposite[0]] = coordinates[collided[0]]  + length_box_x
  
  return(coordinates)