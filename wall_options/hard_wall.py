import numpy as np
from scipy.special import erfinv
import h5py
import arrayfire as af  
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
arrayfire_backend    = params.arrayfire_backend
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

def wall_y(sol):

  y_coordinates  = sol[no_of_particles:2*no_of_particles]

  if(simulation_dimension == 2): 
    velocity_y = sol[3*no_of_particles:4*no_of_particles]

  else:
    velocity_y = sol[4*no_of_particles:5*no_of_particles]
  
  collided_top = np.where(y_coordinates>top_boundary)
  collided_bot = np.where(y_coordinates<bottom_boundary)

  y_coordinates[collided_top[0]] = top_boundary
  velocity_y[collided_top[0]]    = velocity_y[collided_top[0]]*(-1)

  if(simulation_dimension == 2):
    sol[collided_top[0]+3*no_of_particles] = velocity_y[collided_top[0]]
    
  else:
    sol[collided_top[0]+4*no_of_particles] = velocity_y[collided_top[0]]
    
  sol[collided_top[0]+no_of_particles] = y_coordinates[collided_top[0]]

  y_coordinates[collided_bot[0]] = bottom_boundary
  velocity_y[collided_bot[0]]    = velocity_y[collided_bot[0]]*(-1)

  if(simulation_dimension == 2):
    sol[collided_bot[0]+3*no_of_particles] = velocity_y[collided_bot[0]]
    
  else:
    sol[collided_bot[0]+4*no_of_particles] = velocity_y[collided_bot[0]]
    
  sol[collided_bot[0]+no_of_particles] = y_coordinates[collided_bot[0]]
  return(sol)

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