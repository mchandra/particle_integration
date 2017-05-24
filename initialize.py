import numpy as np
from scipy.special import erfinv
import h5py
import params

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

""" 
This file is used in providing the initial velocities and positions to the particles.
The distribution of your choice may be obtained by modifying the options that have been 
provided below. Although the choice of function must be changed from this file, the parameter 
change may also be made at params.py
"""

""" Initializing the positions for the particles """

initial_position_x = left_boundary   + length_box_x * np.random.rand(no_of_particles)
initial_position_y = bottom_boundary + length_box_y * np.random.rand(no_of_particles)
initial_position_z = back_boundary   + length_box_z * np.random.rand(no_of_particles)

""" Initializing velocities to the particles """

# Declaring the random variable which shall be used to sample velocities:
R1 = np.random.rand(no_of_particles)
R2 = np.random.rand(no_of_particles)
R3 = np.random.rand(no_of_particles)

# Sampling velocities corresponding to Maxwell-Boltzmann distribution at T_initial
initial_velocity_x = np.sqrt(2*boltzmann_constant*T_initial/mass_particle)*erfinv(2*R1-1)
initial_velocity_y = np.sqrt(2*boltzmann_constant*T_initial/mass_particle)*erfinv(2*R2-1)
initial_velocity_z = np.sqrt(2*boltzmann_constant*T_initial/mass_particle)*erfinv(2*R3-1)

""" Time parameters for the simulation """

box_crossing_time_scale = (length_box_x/np.max(initial_velocity_x))
final_time              = 5.0  #20 * box_crossing_time_scale
dt                      = 0.001 * box_crossing_time_scale
time                    = np.arange(0, final_time, dt)

""" Writing the data to a file """

if(simulation_dimension == 3):
  h5f = h5py.File('data_files/initial_conditions/initial_data.h5', 'w')
  h5f.create_dataset('time',          data = time)
  h5f.create_dataset('x_coordinates', data = initial_position_x)
  h5f.create_dataset('y_coordinates', data = initial_position_y)
  h5f.create_dataset('z_coordinates', data = initial_position_z)
  h5f.create_dataset('velocity_x',    data = initial_velocity_x)
  h5f.create_dataset('velocity_y',    data = initial_velocity_y)
  h5f.create_dataset('velocity_z',    data = initial_velocity_z)
  h5f.close()

if(simulation_dimension == 2):
  h5f = h5py.File('data_files/initial_conditions/initial_data.h5', 'w')
  h5f.create_dataset('time',          data = time)
  h5f.create_dataset('x_coordinates', data = initial_position_x)
  h5f.create_dataset('y_coordinates', data = initial_position_y)
  h5f.create_dataset('velocity_x',    data = initial_velocity_x)
  h5f.create_dataset('velocity_y',    data = initial_velocity_y)
  h5f.close()
