from initialization import *
import numpy as np
import pylab as pl
import h5py
import scipy.stats as stats

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
    pl.ylim(0,20)
    print ("Time index = ", time_index)
    pl.savefig('temp/%04d'%time_index + '.png')
    pl.clf()
