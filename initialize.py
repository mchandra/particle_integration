import modules
from simulation_parameters import *
from modules import *

""" Initializing the positions for the particles """

initial_position_x = length_box_x * np.random.rand(no_of_particles)
initial_position_y = length_box_y * np.random.rand(no_of_particles)
initial_position_z = length_box_z * np.random.rand(no_of_particles)

""" Initializing velocities to the particles """

# Declaring the random variable which shall be used to sample velocities:
R1 = np.random.rand(no_of_particles)
R2 = np.random.rand(no_of_particles)
R3 = np.random.rand(no_of_particles)

# Sampling velocities corresponding to Maxwell-Boltzmann distribution at T_initial
initial_velocity_x = np.sqrt(2*T_initial)*erfinv(2*R1-1)
initial_velocity_y = np.sqrt(2*T_initial)*erfinv(2*R2-1)
initial_velocity_z = np.sqrt(2*T_initial)*erfinv(2*R3-1)

""" Time parameters for the simulation """

box_crossing_time_scale = (length_box_x/np.max(initial_velocity_x))
final_time              = 60.0 * box_crossing_time_scale
dt                      = 0.005 * box_crossing_time_scale
time                    = np.arange(0, final_time, dt)

""" Combining initial conditions into a single vector """

if(simulation_dimension == 2):
  initial_conditions = np.concatenate([initial_position_x, initial_position_y, \
                                       initial_velocity_y, initial_velocity_y, \
                                      ], axis = 0 \
                                     )

else:
  initial_conditions = np.concatenate([initial_position_x, initial_position_y, initial_position_z, \
                                       initial_velocity_y, initial_velocity_y, initial_velocity_z, \
                                      ], axis = 0 \
                                     )

h5f = h5py.File('data_files/initial_conditions/initial_data.h5', 'w')
h5f.create_dataset('time',                data = time)
h5f.create_dataset('initial_conditions',  data = initial_conditions)
h5f.close()
