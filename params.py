no_of_particles      = 100000
mass_particle        = 1.0
boltzmann_constant   = 1.0
# Change this value to 2 or 3 depending upon whether simulation is in 2D or 3D
simulation_dimension = 2 

# Note that simulation time parameters must be changed from the initialize file

"""Options in development"""

restart_simulation   = "false"
if(restart_simulation  == "true"):
  restart_time_index = 100
  

"""End of options in current developement"""

# If the user requests the written data files will contain the spacial temperature array of each time-step

plot_spatial_temperature_profile = "false"
# Currently the spatial temperature profile is limited to 2D cases

if(plot_spatial_temperature_profile == "true"):
  x_zones = 100
  y_zones = 100

"""Defining the choice for the integrator"""
# The following options are available for integrators:
# Option "verlet"    - Velocity Verlet algorithm
# Option "4th-order" - 4th order symplectic method
# Option "6th-order" - 6th order symplectic method(to be added)
# Option "8th-order" - 8th order symplectic method(to be added)

choice_integrator = "verlet"

""" Defining the collision parameters """
# The following collision kernels are available for implementation
# Option "collisionless"   - Collisionless Model
# Option "hardsphere"      - HS Model
# Option "potential-based" - Potential Based Scattering
# Option "montecarlo"      - MonteCarlo Scattering

collision_operator = "montecarlo"

# scattering_distance is the distance between 2 particles below which scattering kernel is activated
if(collision_operator == "hardsphere"):
  scattering_distance = 0.01

# unlike HS model, the minimum scattering distance can only be altered by shifting the potential function found in collision_operators/potential.py
# default function for potential = potential_magnitude * (-tanh(potential_gradient*distance) + 1)
# here you can change the potential gradient and potential magnitude
elif(collision_operator == "potential-based"):
  potential_steepness     = 300
  potential_amplitude     = 20
  order_finite_difference = 4

# in Monte-Carlo scattering, particles are scattered depending upon the local temperature of the zone they prevail in
# this eliminates the distance check between particles
elif(collision_operator == "montecarlo"):
  x_zones            = 50
  y_zones            = 50
  scattered_fraction = 1.0 # Yet to be implemented. All particles are being scattered currently   

""" Wall and temperature parameters """

# we shall define the different wall conditions that may be implemented:
# Option "periodic" - Periodic B.C's at the walls
# Option "hardwall" - Hardwall B.C's at the walls
# Option "thermal"  - Thermal B.C's at the walls

T_initial        = 1.5

wall_condition_x = "periodic"
wall_condition_y = "periodic"
wall_condition_z = "periodic"

if(wall_condition_x == "thermal"):
  T_left_wall  = 2.0
  T_right_wall = 2.0

if(wall_condition_y == "thermal"):
  T_top_wall = 2.0
  T_bot_wall = 2.0

if(wall_condition_z == "thermal"):
  T_front_wall = 2.0
  T_back_wall  = 2.0

 
""" Length Parameters of Simulation Domain """

left_boundary    = 0.
right_boundary   = 1.
length_box_x     = right_boundary - left_boundary

bottom_boundary  = 0.
top_boundary     = 1.
length_box_y     = top_boundary   - bottom_boundary

back_boundary    = 0.
front_boundary   = 1.
length_box_z     = front_boundary - back_boundary