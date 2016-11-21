from initialization import *
import numpy as np
import pylab as pl
import h5py
import scipy.stats as stats

"""Plotting of the temperature distribution in each zome"""

temperature=np.zeros((x.size-1,time.size),dtype=np.float)

for time_index,t0 in enumerate(time):
    if(time_index==time.size-1):
        break
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'r')
    solution_all_current = h5f['solution_all/solution_dataset_'+str(time_index)][:]
    h5f.close()
    print("Temperature Computation for time_index : "+time_index)
    for k in range(0,x.size-1)
        for p in range(0,no_of_particles):
            for i in range(0,x.size-1):
                if((solution_all_current[p]>x[i])and(solution_all_current[p]<x[i+1])):
                    x_zone = i
            for j in range(0,y.size-1):
                if((solution_all_current[p+no_of_particles]>y[j])and(solution_all_current[p+no_of_particles]<y[j+1])):
                    y_zone = j
            
            if(x_zone==k):
                temperature[k,time_index]=solution_all_current[p+2*no_of_particles]**2+solution_all_current[p+3*no_of_particles]**2

"""Creating Density Plot Movie"""
for time_index,t0 in enumerate(time[::10]):          
    pl.xlabel('$x$')
    pl.ylabel('$\mathrm{Temperature}$')
    pl.plot(temperature[:,time_index])
    print ("Time index = ", time_index)
    pl.savefig('temperature/%04d'%time_index + '.png')
    pl.clf()
