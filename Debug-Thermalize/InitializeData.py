import numpy as np
import h5py
from SimulationData import *

""" Initial conditions """


initial_conditions_position_x = left_boundary + length_of_box_x*np.random.rand(no_of_particles)
initial_conditions_position_y = bottom_boundary + length_of_box_y*np.random.rand(no_of_particles)
initial_conditions_position_z = back_boundary + length_of_box_z*np.random.rand(no_of_particles)
initial_conditions_velocity_x=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_y=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_z=np.zeros(no_of_particles,dtype=np.float)

for i in range(0,no_of_particles,2):

    x_1=np.random.rand(1)
    x_2=np.random.rand(1)
    y_1=const*np.sqrt(-2*np.log(x_1))*np.cos(2*np.pi*x_2)
    y_2=const*np.sqrt(-2*np.log(x_1))*np.sin(2*np.pi*x_2)
    x_3=np.random.rand(1)
    x_4=np.random.rand(1)
    y_3=const*np.sqrt(-2*np.log(x_3))*np.cos(2*np.pi*x_4)
    y_4=const*np.sqrt(-2*np.log(x_3))*np.sin(2*np.pi*x_4)
    x_5=np.random.rand(1)
    x_6=np.random.rand(1)
    y_5=const*np.sqrt(-2*np.log(x_5))*np.cos(2*np.pi*x_6)
    y_6=const*np.sqrt(-2*np.log(x_5))*np.sin(2*np.pi*x_6)
    initial_conditions_velocity_x[i]=y_1
    initial_conditions_velocity_y[i]=y_2
    initial_conditions_velocity_z[i]=y_3
    initial_conditions_velocity_x[i+1]=y_4
    initial_conditions_velocity_y[i+1]=y_5
    initial_conditions_velocity_z[i+1]=y_6


""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                                    initial_conditions_position_z, initial_conditions_velocity_x,\
                                    initial_conditions_velocity_y,initial_conditions_velocity_z],axis = 0)

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)
final_time            = 50 * box_crossing_time_scale
dt   = 0.05 * box_crossing_time_scale
time = np.arange(0, final_time, dt)

"""Creating the intial conditions file"""

h5f = h5py.File('Initial.h5', 'w')
h5f.create_dataset('time', data=time)
h5f.create_dataset('initial_conditions', data=initial_conditions)
h5f.close() 


