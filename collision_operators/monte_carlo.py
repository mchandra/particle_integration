import numpy as np
from scipy.special import erfinv
import h5py
import arrayfire as af  
import params

"""
This script incorporates the montecarlo method of scattering.
By this method we divide our region of interest into zones and 
scatter the particles using an MB distribution, corresponding to 
the temperature of that zone
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

def collision_operator(sol):
  
  x_zones_particle   = (x_zones/length_box_x) * sol[0:no_of_particles]
  y_zones_particle   = (y_zones/length_box_y) * sol[no_of_particles:2*no_of_particles]
  x_zones_particle   = x_zones_particle.astype(int)
  y_zones_particle   = y_zones_particle.astype(int)
  zone               = x_zones*y_zones_particle + x_zones_particle
  zonecount          = np.bincount(zone)
  
  temp = np.zeros(zonecount.size)

  for i in range(x_zones*y_zones):
    indices = np.where(zone == i)[0]
    temp[i] = 0.5*np.sum(sol[indices+2*no_of_particles]**2 + sol[indices+3*no_of_particles]**2)
  
  temp=temp/zonecount
    
  for i in range(x_zones*y_zones):
    indices = np.where(zone == i)[0]
    x1 = np.random.rand(zonecount[i])
    x2 = np.random.rand(zonecount[i])
    sol[indices+2*no_of_particles] = np.sqrt(2*temp[i])*erfinv(2*x1-1) 
    sol[indices+3*no_of_particles] = np.sqrt(2*temp[i])*erfinv(2*x2-1)

  return(sol)