import numpy as np
from scipy.special import erfinv
import h5py
import params

""" 
This file contains 3 functions, which define hardwall B.C's in 3 directions
Depending upon the choice of the user, hardwall boundary conditions may be set
to either of the x,y and z directions.

A hardwall B.C means that a particle that encounters such a boundary will reflect
back with its component of velocity perpendicular to the wall being the same in
magnitude, but opposite in direction
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

def wall_x(x_coords, vel_x, vel_y, vel_z):

  collided_right = np.where(x_coords > right_boundary)
  collided_left  = np.where(x_coords < left_boundary)

  x_coords[collided_left[0]] = left_boundary
  vel_x[collided_left[0]]    = vel_x[collided_left[0]]*(-1)    

  x_coords[collided_right[0]] = right_boundary
  vel_x[collided_right[0]]    = vel_x[collided_right[0]]*(-1)    

  return(x_coords, vel_x, vel_y, vel_z)

def wall_y(y_coords, vel_x, vel_y, vel_z):

  collided_top = np.where(y_coords > top_boundary)
  collided_bot = np.where(y_coords < bot_boundary)

  y_coords[collided_bot[0]] = bot_boundary
  vel_y[collided_bot[0]]    = vel_y[collided_bot[0]]*(-1)    

  y_coords[collided_top[0]] = top_boundary
  vel_y[collided_top[0]]    = vel_y[collided_top[0]]*(-1)    

  return(y_coords, vel_x, vel_y, vel_z)

def wall_z(z_coords, vel_x, vel_y, vel_z):

  collided_front = np.where(z_coords > front_boundary)
  collided_back  = np.where(z_coords < back_boundary)

  z_coords[collided_back[0]] = back_boundary
  vel_z[collided_back[0]]    = vel_z[collided_back[0]]*(-1)    

  z_coords[collided_front[0]] = front_boundary
  vel_z[collided_front[0]]    = vel_z[collided_front[0]]*(-1)    

  return(z_coords, vel_x, vel_y, vel_z)