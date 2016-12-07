import numpy as np
import h5py


no_of_particles=1000000
length_of_box_x=1
dx = (length_of_box_x/32)
x = np.arange(0,1,dx)
x= np.concatenate((x,[1]),axis = 0)

h5f = h5py.File('initial_conditions.h5', 'r')
initial_conditions = h5f['initial_conditions_dataset'][:]
h5f.close()
v_x=np.zeros(no_of_particles,dtype=np.float)
v_x=initial_conditions[3*no_of_particles:4*no_of_particles]

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.amax(v_x)
final_time            = 4 * box_crossing_time_scale
dt   = 0.005 * box_crossing_time_scale
time = np.arange(0, final_time, dt)

""" Computing density plots """

amp = np.zeros(time.size,dtype=np.float)

for time_index,t0 in enumerate(time):
    if(time_index==time.size-1):
        break
    dumpfile   = 'solution_all/solution_' + str(time_index) + '.h5'
    datafile   = h5py.File(dumpfile, "r")
    solution   = datafile['solution_all']['solution_dataset_' + str(time_index)]
    n, bins = np.histogram(solution[:no_of_particles], bins=x)
    amp[time_index]=np.amax(n)

amp=amp/no_of_particles
amp=32*amp
amp=amp-1    
h5f = h5py.File('plot.h5', 'w')
h5f.create_dataset('time', data=time)
h5f.create_dataset('amp', data=amp)
h5f.close()

