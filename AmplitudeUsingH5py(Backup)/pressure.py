import numpy as np
import pylab as pl
import h5py
import scipy.stats as stats

no_of_particles = 100000
length_of_box_x = 1
h5f = h5py.File('initial_conditions.h5', 'r')
initial_conditions = h5f['initial_conditions_dataset'][:]
h5f.close()
print(initial_conditions,"	",initial_conditions.size)    
v_x=np.zeros(no_of_particles,dtype=np.float)
v_x=initial_conditions[2*no_of_particles:3*no_of_particles]

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(v_x)
final_time            = 5 * box_crossing_time_scale
dt   = 0.01 * box_crossing_time_scale
time = np.arange(0, final_time, dt)   
x=np.linspace(0,1,33)

for time_index,t0 in enumerate(time):
    if(time_index==time.size-1):
        break
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'r')
    solution_all_current = h5f['solution_all/solution_dataset_'+str(time_index)][:]
    h5f.close()
    pressure=0
    for i in range(no_of_particles):
        pressure=pressure+solution_all_current[i+2*no_of_particles]**2+solution_all_current[i+3*no_of_particles]**2
    pressure=pressure/no_of_particles
    print("Pressure = ",pressure)
    
