from solver import *
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



for time_index,t0 in enumerate(time[::10]):
    for i in range(0,no_of_particles):
        if(time_index==time.size-1):
            break
        h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'r')
        s = h5f['solution_all/solution_dataset_'+str(time_index)][:]
        h5f.close()
        pl.plot(s[i],s[no_of_particles + i], 'o',color='blue', markersize=10, alpha = 0.4)
        pl.axhline(y=-1)
        pl.axvline(x=-1)
        pl.axhline(y=1)
        pl.axvline(x=1)
        pl.xlim([left_boundary, right_boundary])
        pl.ylim([bottom_boundary, top_boundary])
        pl.title('$\mathrm{Time}$ = ' + str(time[time_index]) )
        pl.xlabel('$x$')
        pl.ylabel('$y$')
    print ("Time index = ", time_index)
    pl.savefig('images/point_mass' + '%04d'%time_index + '.png')
    pl.clf() 
