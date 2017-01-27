no_of_particles      = 100000
mass_particle        = 1.0
boltzmann_constant   = 1.0
simulation_dimension = 2 # Change this value to 2 or 3 depending upon the need

"""Defining the choice for the integrator"""
# The following options are available for integrators:
# Option 0 - Velocity Verlet algorithm
# Option 1 - 4th order symplectic method
choice_integrator = 0

""" Defining the collision parameters """
# The following collision kernels are available for implementation
# Option 0 - Collisionless Model
# Option 1 - HS Model
# Option 2 - Potential Based Scattering
# Option 3 - MonteCarlo Scattering

collision_operator = 0

# scattering_distance is the distance between 2 particles below which scattering kernel is activated
if(collision_operator == 1):
  scattering_distance = 0.01

# unlike HS model, the minimum scattering distance can only be altered by shifting the potential function found in collision_operators/potential.py
# default function for potential = potential_magnitude * (-tanh(potential_gradient*distance) + 1)
# here you can change the potential gradient and potential magnitude
elif(collision_operator == 2):
  potential_steepness     = 300
  potential_amplitude     = 20
  order_finite_difference = 4

# in Monte-Carlo scattering, particles are scattered depending upon the local temperature of the zone they prevail in
# this eliminates the distance check between particles
elif(collision_operator == 3):
  x_zones            = 10
  y_zones            = 10
  scattered_fraction = 1.0   

""" Wall and temperature parameters """

# we shall define the different wall conditions that may be implemented:
# Option 0 - Periodic B.C's at the walls
# Option 1 - Hardwall B.C's at the walls
# Option 2 - Thermal B.C's at the walls

T_initial    = 1.5

wall_condition_x = 2
wall_condition_y = 0
wall_condition_z = 0

if(wall_condition_x == 2):
  T_left_wall  = 2.0
  T_right_wall = 2.0

if(wall_condition_y == 2):
  T_top_wall = 2.0
  T_bot_wall = 2.0

if(wall_condition_z == 2):
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