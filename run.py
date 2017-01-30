import numpy as np
from scipy.special import erfinv
import h5py
import arrayfire as af  
import params

"""Here we shall assign values as set in params"""

no_of_particles      = params.no_of_particles
simulation_dimension = params.simulation_dimension
restart_simulation   = params.restart_simulation
arrayfire_backend    = params.arrayfire_backend
choice_integrator    = params.choice_integrator
collision_operator   = params.collision_operator

if(restart_simulation == "true"):
  restart_time_index = params.restart_time_index

plot_spatial_temperature_profile = params.plot_spatial_temperature_profile

if(plot_spatial_temperature_profile == "true"):
  x_zones = params.x_zones
  y_zones = params.y_zones

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

if(simulation_dimension == 2):
  sol = np.zeros(4*no_of_particles,dtype=np.float)
  old = np.zeros(4*no_of_particles,dtype=np.float)

if(simulation_dimension == 3):
  sol = np.zeros(6*no_of_particles,dtype=np.float)
  old = np.zeros(6*no_of_particles,dtype=np.float)

h5f                = h5py.File('data_files/initial_conditions/initial_data.h5', 'r')
initial_conditions = h5f['initial_conditions'][:]
time               = h5f['time'][:]
h5f.close()

if(restart_simulation == "true"):
  h5f                = h5py.File('data_files/timestepped_data/solution_'+str(restart_time_index)+'.h5', 'r')
  initial_conditions = h5f['sol'][:]
  time               = h5f['time'][:]
  h5f.close()

dt = time[1]

"""Declaring data variables which shall be used in post-processing"""

momentum_x     = np.zeros(time.size)
momentum_y     = np.zeros(time.size)
kinetic_energy = np.zeros(time.size)
pressure       = np.zeros(time.size)
heatflux_x     = np.zeros(time.size)
heatflux_y     = np.zeros(time.size)

if(simulation_dimension == 3):
  momentum_z = np.zeros(time.size)
  heatflux_z = np.zeros(time.size)

if(collision_operator == "potential-based"):
  potential_energy = np.zeros(time.size)

"""Choice for integrators"""
if(choice_integrator == "verlet"):
  from integrators.verlet import integrator

# Not fully stable. Few bugs may be present
elif(choice_integrator == "4th-order"):
  from integrators.fourth_order_symplectic import integrator

"""Setting the wall options"""

if(wall_condition_x == "thermal"):
  from wall_options.thermal  import wall_x
elif(wall_condition_x == "hardwall"):
  from wall_options.hard_wall import wall_x
elif(wall_condition_x == "periodic"):
  from wall_options.periodic import wall_x

if(wall_condition_y == "thermal"):
  from wall_options.thermal import wall_y
elif(wall_condition_y == "hardwall"):
  from wall_options.hard_wall import wall_y
elif(wall_condition_y == "periodic"):
  from wall_options.periodic import wall_y

if(simulation_dimension == 3):
  if(wall_condition_z == "thermal"):
    from wall_options.thermal import wall_z
  elif(wall_condition_z == "hardwall"):
    from wall_options.hard_wall import wall_z
  elif(wall_condition_z == "periodic"):
    from wall_options.periodic import wall_z

"""Collision Options"""

if(collision_operator == "hardsphere"):
  from collision_operators.hard_sphere import collision_operator

if(collision_operator == "montecarlo"):
  from collision_operators.monte_carlo import collision_operator

# Now we shall proceed to evolve the system with time:

for time_index,t0 in enumerate(time):
  
  if(restart_simulation == "true"):
    time_index = time_index + restart_time_index

  print("Computing For Time Index = ",time_index)
  
  
  if(time_index == time.size-1):
    break
  
  if(time_index==0):
    initial_conditions = initial_conditions
  else:
    initial_conditions = old
  
  sol = integrator(initial_conditions,dt)
  sol = wall_x(sol)
  sol = wall_y(sol)
  
  if(simulation_dimension == 3):
    sol = wall_z(sol)
  
  if(collision_operator == "hardwall" or collision_operator == "montecarlo"):
    sol = collision_operator(sol)
  
  old = sol

  """Declaring variables used in calculation for post-processor"""

  if(simulation_dimension == 2):
    x_coordinates = sol[0:no_of_particles]
    y_coordinates = sol[no_of_particles:2*no_of_particles] 
    velocity_x    = sol[2*no_of_particles:3*no_of_particles]
    velocity_y    = sol[3*no_of_particles:4*no_of_particles]

  if(simulation_dimension == 3):
    x_coordinates = sol[0:no_of_particles]
    y_coordinates = sol[no_of_particles:2*no_of_particles] 
    z_coordinates = sol[2*no_of_particles:3*no_of_particles]
    velocity_x    = sol[3*no_of_particles:4*no_of_particles]
    velocity_y    = sol[4*no_of_particles:5*no_of_particles]
    velocity_z    = sol[5*no_of_particles:6*no_of_particles]

  if(collision_operator == "montecarlo"):
    from collision_operators.monte_carlo import collision_operator
    sol = collision_operator(sol)

  if(plot_spatial_temperature_profile == "true" and time_index%100 == 0):
    particle_xzone   = (x_zones/length_box_x) * sol[0:no_of_particles]
    particle_yzone   = (y_zones/length_box_y) * sol[no_of_particles:2*no_of_particles]
    particle_xzone   = particle_xzone.astype(int)
    particle_yzone   = particle_yzone.astype(int)
    particle_zone    = x_zones * particle_yzone + particle_xzone
    zonecount        = np.bincount(particle_zone)
    
    temperature_spatial = np.zeros(zonecount.size)
    

    for i in range(x_zones*y_zones):
      indices = np.where(particle_zone == i)[0]
      temperature_spatial[i] = 0.5*np.sum(sol[indices+2*no_of_particles]**2 + sol[indices+3*no_of_particles]**2)
    
    temperature_spatial = temperature_spatial/zonecount
    temperature_spatial = temperature_spatial.reshape(x_zones,y_zones)
  
  """Calculation of the functions which will be used to post-process the results of the simulation run"""

  if(simulation_dimension == 2):
    momentum_x[time_index]     = mass_particle * np.sum(velocity_x)
    momentum_y[time_index]     = mass_particle * np.sum(velocity_y)
    kinetic_energy[time_index] = 0.5*mass_particle*np.sum(velocity_x**2 + velocity_y**2)
    pressure[time_index]       = np.sum(velocity_x**2 + velocity_y**2)/no_of_particles
    heatflux_x[time_index]     = np.sum(velocity_x*(velocity_x**2 + velocity_y**2))/no_of_particles
    heatflux_y[time_index]     = np.sum(velocity_y*(velocity_x**2 + velocity_y**2))/no_of_particles

  if(simulation_dimension == 3):
    momentum_x[time_index]     = mass_particle * np.sum(velocity_x)
    momentum_y[time_index]     = mass_particle * np.sum(velocity_y)
    momentum_z[time_index]     = mass_particle * np.sum(velocity_z)

    kinetic_energy[time_index] = 0.5*mass_particle*np.sum(velocity_x**2 + velocity_y**2 + velocity_z**2)
    
    pressure[time_index]       = np.sum(velocity_x**2 + velocity_y**2 + velocity_z**2)/no_of_particles
    
    heatflux_x[time_index]     = np.sum(velocity_x*(velocity_x**2 + velocity_y**2 + velocity_z**2))/no_of_particles
    heatflux_y[time_index]     = np.sum(velocity_y*(velocity_x**2 + velocity_y**2 + velocity_z**2))/no_of_particles
    heatflux_z[time_index]     = np.sum(velocity_z*(velocity_x**2 + velocity_y**2 + velocity_z**2))/no_of_particles

  if(collision_operator == "potential-based"):
    from collision_operators.potential import calculate_potential_energy
    potential_energy = calculate_potential_energy(sol)

  print(pressure[time_index])
  
  # Writing the data to file every 1000 time steps
  # This data will then be post-processed to generate results
  
  if((time_index%100)==0):
    h5f = h5py.File('data_files/timestepped_data/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('sol',                data = sol)
    h5f.create_dataset('momentum_x',         data = momentum_x)
    h5f.create_dataset('momentum_y',         data = momentum_y)
    h5f.create_dataset('heatflux_x',         data = heatflux_x)
    h5f.create_dataset('heatflux_y',         data = heatflux_y)
    if(simulation_dimension == 3):
      h5f.create_dataset('heatflux_z',       data = heatflux_z)
      h5f.create_dataset('momentum_z',       data = momentum_z)
    h5f.create_dataset('kinetic_energy',     data = kinetic_energy)
    h5f.create_dataset('pressure',           data = pressure)
    if(collision_operator == "potential-based"):
      h5f.create_dataset('potential_energy', data = potential_energy)
    if(plot_spatial_temperature_profile == "true"):
      h5f.create_dataset('temperature_spatial', data = temperature_spatial)
    h5f.create_dataset('time',               data = time)
    h5f.close()