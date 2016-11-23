
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



for time_index,t0 in enumerate(time):
    if(time_index==time.size-1):
        break
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'r')
    solution_all_current = h5f['solution_all/solution_dataset_'+str(time_index)][:]
    h5f.close()
    print ("Time index = ", time_index)
    for p in range(no_of_particles):
        pl.title('Numerical Density along center line with time')
        pl.xlabel('$x$')
        pl.ylabel('$\mathrm{velocity_x}$')
        pl.plot(solution_all_current[p],solution_all_current[p+(2*no_of_particles)], '.',color='blue', markersize=6, alpha = 0.4)
        
    pl.savefig('post/phase/%04d'%time_index + '.png')
    pl.clf()
       
 
