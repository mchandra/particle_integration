import numpy as np
from scipy.special import erfinv
import h5py 
import params

""" 
This file contains 3 functions, which define thermal B.C's in 3 directions
Depending upon the choice of the user, thermal boundary conditions may be set
to either of the x,y and z directions.

A thermal B.C means that a particle that encounters such a boundary will reflect
back with its component of vel perpendicular to the wall, away from it with 
its magnitude taking a value corresponding to the temperature of the wall
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

if(simulation_dimension == 3):

  def wall_x(x_coords, vel_x, vel_y, vel_z):

    collided_right = np.where(x_coords>right_boundary)
    collided_left  = np.where(x_coords<left_boundary)

    # Random variables used to assign the new velocities to the particles:
    
    R1 = np.random.rand(collided_right[0].size)
    R2 = np.random.rand(collided_right[0].size)
    R3 = np.random.rand(collided_right[0].size)
    
    x_coords[collided_right[0]] = right_boundary - 1e-12

    vel_x[collided_right[0]] = np.sqrt(-2*T_right_wall*(boltzmann_constant/mass_particle)*np.log(R1))*(-1)    
    vel_y[collided_right[0]] = np.sqrt(2*T_right_wall*(boltzmann_constant/mass_particle))*erfinv(2*R2-1)
    vel_z[collided_right[0]] = np.sqrt(2*T_right_wall*(boltzmann_constant/mass_particle))*erfinv(2*R3-1)

    R1 = np.random.rand(collided_left[0].size)
    R2 = np.random.rand(collided_left[0].size)
    R3 = np.random.rand(collided_left[0].size)

    x_coords[collided_left[0]] = left_boundary - 1e-12

    vel_x[collided_left[0]] = np.sqrt(-2*T_left_wall*(boltzmann_constant/mass_particle)*np.log(R1))   
    vel_y[collided_left[0]] = np.sqrt(2*T_left_wall*(boltzmann_constant/mass_particle))*erfinv(2*R2-1)
    vel_z[collided_left[0]] = np.sqrt(2*T_left_wall*(boltzmann_constant/mass_particle))*erfinv(2*R3-1)

    return(x_coords, vel_x, vel_y, vel_z)

  def wall_y(y_coords, vel_x, vel_y, vel_z):

    y_coords  = sol[no_of_particles:2*no_of_particles]

    collided_top = np.where(y_coords>top_boundary)
    collided_bot = np.where(y_coords<bottom_boundary)

    # Random variables used to assign the new velocities to the particles:
    R1 = np.random.rand(collided_top[0].size)
    R2 = np.random.rand(collided_top[0].size)
    R3 = np.random.rand(collided_top[0].size)
    
    y_coords[collided_top[0]] = top_boundary - 1e-12

    vel_x[collided_top[0]] = np.sqrt(2*T_top_wall*(boltzmann_constant/mass_particle))*erfinv(2*R1-1)    
    vel_y[collided_top[0]] = np.sqrt(-2*T_top_wall*(boltzmann_constant/mass_particle)*np.log(R2))*(-1)
    vel_z[collided_top[0]] = np.sqrt(2*T_top_wall*(boltzmann_constant/mass_particle))*erfinv(2*R3-1)

    R1 = np.random.rand(collided_bot[0].size)
    R2 = np.random.rand(collided_bot[0].size)
    R3 = np.random.rand(collided_bot[0].size)

    y_coords[collided_bot[0]] = bottom_boundary - 1e-12

    vel_x[collided_bot[0]] = np.sqrt(2*T_bot_wall*(boltzmann_constant/mass_particle))*erfinv(2*R1-1)    
    vel_y[collided_bot[0]] = np.sqrt(-2*T_bot_wall*(boltzmann_constant/mass_particle)*np.log(R2))
    vel_z[collided_bot[0]] = np.sqrt(2*T_bot_wall*(boltzmann_constant/mass_particle))*erfinv(2*R3-1)

    return(y_coords, vel_x, vel_y, vel_z)

if(simulation_dimension == 2):

  def wall_x(x_coords, vel_x, vel_y):

    collided_right = np.where(x_coords>right_boundary)
    collided_left  = np.where(x_coords<left_boundary)

    # Random variables used to assign the new velocities to the particles:
    
    R1 = np.random.rand(collided_right[0].size)
    R2 = np.random.rand(collided_right[0].size)
    
    x_coords[collided_right[0]] = right_boundary - 1e-12

    vel_x[collided_right[0]] = np.sqrt(-2*T_right_wall*(boltzmann_constant/mass_particle)*np.log(R1))*(-1)    
    vel_y[collided_right[0]] = np.sqrt(2*T_right_wall*(boltzmann_constant/mass_particle))*erfinv(2*R2-1)

    R1 = np.random.rand(collided_left[0].size)
    R2 = np.random.rand(collided_left[0].size)

    x_coords[collided_left[0]] = left_boundary - 1e-12

    vel_x[collided_left[0]] = np.sqrt(-2*T_left_wall*(boltzmann_constant/mass_particle)*np.log(R1))   
    vel_y[collided_left[0]] = np.sqrt(2*T_left_wall*(boltzmann_constant/mass_particle))*erfinv(2*R2-1)

    return(x_coords, vel_x, vel_y)

  def wall_y(y_coords, vel_x, vel_y):

    collided_top = np.where(y_coords>top_boundary)
    collided_bot = np.where(y_coords<bottom_boundary)

    # Random variables used to assign the new velocities to the particles:
    R1 = np.random.rand(collided_top[0].size)
    R2 = np.random.rand(collided_top[0].size)
    
    y_coords[collided_top[0]] = top_boundary - 1e-12

    vel_x[collided_top[0]] = np.sqrt(2*T_top_wall*(boltzmann_constant/mass_particle))*erfinv(2*R1-1)    
    vel_y[collided_top[0]] = np.sqrt(-2*T_top_wall*(boltzmann_constant/mass_particle)*np.log(R2))*(-1)

    R1 = np.random.rand(collided_bot[0].size)
    R2 = np.random.rand(collided_bot[0].size)

    y_coords[collided_bot[0]] = bottom_boundary - 1e-12

    vel_x[collided_bot[0]] = np.sqrt(2*T_bot_wall*(boltzmann_constant/mass_particle))*erfinv(2*R1-1)    
    vel_y[collided_bot[0]] = np.sqrt(-2*T_bot_wall*(boltzmann_constant/mass_particle)*np.log(R2))

    return(y_coords, vel_x, vel_y)

def wall_z(z_coords, vel_x, vel_y, vel_z):

  collided_front = np.where(z_coords>front_boundary)
  collided_back  = np.where(z_coords<back_boundary)

  # Random variables used to assign the new velocities to the particles:
  R1 = np.random.rand(collided_front[0].size)
  R2 = np.random.rand(collided_front[0].size)
  R3 = np.random.rand(collided_front[0].size)

  z_coords[collided_front[0]] = front_boundary - 1e-12

  vel_x[collided_front[0]] = np.sqrt(2*T_front_wall*(boltzmann_constant/mass_particle))*erfinv(2*R1-1)    
  vel_y[collided_front[0]] = np.sqrt(2*T_front_wall*(boltzmann_constant/mass_particle))*erfinv(2*R2-1)
  vel_z[collided_front[0]] = np.sqrt(-2*T_front_wall*(boltzmann_constant/mass_particle)*np.log(R3))*(-1)

  R1 = np.random.rand(collided_back[0].size)
  R2 = np.random.rand(collided_back[0].size)
  R3 = np.random.rand(collided_back[0].size)

  z_coords[collided_back[0]] = back_boundary - 1e-12

  vel_x[collided_back[0]] = np.sqrt(2*T_back_wall*(boltzmann_constant/mass_particle))*erfinv(2*R1-1)    
  vel_y[collided_back[0]] = np.sqrt(2*T_back_wall*(boltzmann_constant/mass_particle))*erfinv(2*R2-1)
  vel_z[collided_back[0]] = np.sqrt(-2*T_back_wall*(boltzmann_constant/mass_particle)*np.log(R3))
  
  return(z_coords, vel_x, vel_y, vel_z)