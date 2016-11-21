
from initialization import *
import numpy as np
import pylab as pl
import h5py


""" Set plot parameters to make beautiful plots """
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'

pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.major.pad']  = 8
pl.rcParams['xtick.minor.pad']  = 8
pl.rcParams['xtick.color']      = 'k'
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'

pl.rcParams['ytick.major.size'] = 8
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.major.pad']  = 8
pl.rcParams['ytick.minor.pad']  = 8
pl.rcParams['ytick.color']      = 'k'
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'




""" Initializing density matrix """

density = np.zeros((x.size-1,y.size-1,time.size),dtype=np.float)




""" Computing density plots """

i=0
j=0
x_zone = 0
y_zone = 0


for time_index,t0 in enumerate(time):
    if(time_index==time.size-1):
        break
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'r')
    solution_all_current = h5f['solution_all/solution_dataset_'+str(time_index)][:]
    h5f.close()
    print("Density Computation for time_index : ",time_index)
    for p in range(0,no_of_particles):
        for i in range(0,x.size-1):
            if((solution_all_current[p]>x[i])and(solution_all_current[p]<x[i+1])):
                x_zone = i
        for j in range(0,y.size-1):
            if((solution_all_current[p+no_of_particles]>y[j])and(solution_all_current[p+no_of_particles]<y[j+1])):
                y_zone = j
        
        density[x_zone,y_zone,time_index] = density[x_zone,y_zone,time_index] +1
       



""" normalizing the density """
for time_index,t0 in enumerate(time):
    density[:,:,time_index] = (density[:,:,time_index]/(no_of_particles/(x_divisions*y_divisions)))
   
amp = np.zeros(time.size)
for time_index,t0 in enumerate(time):
    amp[time_index]=np.amax(density[:,:,time_index])

 

"""Creating Density Plot Movie"""
         
pl.title('Amplitude of Numerical Density')
pl.xlabel('$x$')
pl.ylabel('$\mathrm{Density}$')
pl.plot(amp)
pl.ylim(1,1.6)
pl.savefig('post/amplitude.png')
pl.clf()

    

"""Creating Density Plot Movie"""
for time_index,t0 in enumerate(time[::10]):          
    pl.title('Numerical Density along center line with time')
    pl.xlabel('$x$')
    pl.ylabel('$\mathrm{Density}$')
    pl.ylim(0.4,1.6)
    pl.plot(density[:,int(y_divisions/2),time_index])
    print ("Time index = ", time_index)
    pl.savefig('post/%04d'%time_index + '.png')
    pl.clf()
