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
x=np.zeros(33,dtype=np.float)

"""Plotting of the temperature distribution in each zome"""

for time_index,t0 in enumerate(time):
    if(time_index==time.size-1):
        break
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'r')
    solution_all_current = h5f['solution_all/solution_dataset_'+str(time_index)][:]
    h5f.close()
    print("Temperature Computation for time_index : ",time_index)
    temperature=np.zeros(x.size-1,dtype=np.float)
    count=np.zeros(x.size-1,dtype=np.float)
    for p in range(0,no_of_particles):
     
        for i in range(0,x.size-1):
            if((solution_all_current[p]>x[i])and(solution_all_current[p]<x[i+1])):
                count[i]=count[i]+1
                temperature[i]=temperature[i]+solution_all_current[p+2*no_of_particles]**2+solution_all_current[p+3*no_of_particles]**2
                     
    for i in range(0,x.size-1):
        temperature[i]=temperature[i]/count[i]    
          
    pl.xlabel('$x$')
    pl.ylabel('$\mathrm{Temperature}$')
    pl.plot(x[0:(x.size-1)],temperature)
    pl.savefig('temp/%04d'%time_index + '.png')
    pl.clf()
