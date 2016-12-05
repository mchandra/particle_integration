
from initialization import *
import numpy as np
import pylab as pl
import h5py
no_of_particles=100000

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


""" Initializing density matrix """

density = np.zeros((x.size-1,y.size-1,time.size),dtype=np.float)




""" Computing density plots """



for time_index,t0 in enumerate(time):
    if(time_index==time.size-1):
        break
    dumpfile   = 'solution_all/solution_' + str(time_index) + '.h5'
    datafile   = h5py.File(dumpfile, "r")
    solution   = datafile['solution_all']['solution_dataset_' + str(time_index)]    


    n, bins, patches = pl.hist(solution[:no_of_particles], bins=x, histtype='step')
    pl.clf()
    pl.plot(x[:-1], n/(no_of_particles/(x_divisions*y_divisions)), 'o-')
    pl.xlim([-0.1, 1.1])
    pl.ylim([0.4,1.6])
    pl.xlabel('$x$')
    pl.ylabel('$\mathrm{Number\;of\;Particles}$')
    pl.savefig('post/%04d'%time_index + '.png')
    pl.clf()    


""" normalizing the density 
for time_index,t0 in enumerate(time):
    density[:,:,time_index] = (density[:,:,time_index])"""
   
amp = np.zeros(time.size)
for time_index,t0 in enumerate(time):
    amp[time_index]=np.amax(density[:,:,time_index])

         
pl.title('Amplitude of Numerical Density')
pl.xlabel('$x$')
pl.ylabel('$\mathrm{Density}$')
pl.plot(amp)
pl.ylim(1,1.6)
pl.savefig('amplitude.png')

    

"""Creating Density Plot Movie
for time_index,t0 in enumerate(time[::10]):          
    pl.title('Numerical Density along center line with time')
    pl.xlabel('$x$')
    pl.ylabel('$\mathrm{Density}$')
    pl.ylim(0.4,1.6)
    pl.plot(density[:,int(y_divisions/2),time_index])
    print ("Time index = ", time_index)"""

