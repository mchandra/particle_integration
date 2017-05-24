import numpy as np
from scipy.special import erfinv
import h5py
import params

"""
This model uses the HS model for scattering the particles.
By this model the particles are effectively scattered in the 
same way as to how one models collisions amongst billiard 
balls. By this model we shall condsider only 2-body collisions.
Multi-body collisions are treated as 2 body collisions, with the 
choice of the 2 bodies being random.
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

#Here we complete import of all the variable from the parameters file'

def collision_operator(x_coords, y_coords, z_coords, vel_x, vel_y, vel_z):
  sigma = np.pi * d**2
  # In the above equation d is the value of the diameter of the particles
  make_square_matrix = np.ones(no_of_particles)
  vel_x              = vel_x*make_square_matrix
  vel_x              = vel_x - vel_x.transpose()
  vel_y              = vel_y*make_square_matrix
  vel_y              = vel_y - vel_y.transpose()
  vel_z              = vel_z*make_square_matrix
  vel_z              = vel_z - vel_z.transpose()
  g                  = np.sqrt(vel_x**2 + vel_y**2 + vel_z**2)
  probability        = sigma*g*dt/(
