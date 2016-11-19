from solver import *
from initialization import *
from solver import *
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

density = np.zeros((x.size-1,y.size-1),dtype=np.float)




""" Computing density plots """

i=0
j=0
x_zone = 0
y_zone = 0



h5f = h5py.File('initial_conditions.h5', 'r')
solution_all_current = h5f['initial_conditions_dataset'][:]
h5f.close()
for p in range(0,no_of_particles):
    for i in range(0,x.size-1):
        if((solution_all_current[p]>x[i])and(solution_all_current[p]<x[i+1])):
            x_zone = i
    for j in range(0,y.size-1):
        if((solution_all_current[p+no_of_particles]>y[j])and(solution_all_current[p+no_of_particles]<y[j+1])):
            y_zone = j
    
    density[x_zone,y_zone] = density[x_zone,y_zone] +1


""" normalizing the density """

density[:,:] = density[:,:]/(no_of_particles/(x_divisions*y_divisions))

"""Creating Density Plot Movie"""
        
pl.title('Density Along Centerline')
pl.xlabel('$x$')
pl.ylabel('$\mathrm{Density}$')
pl.ylim(0,2)
pl.plot(density[:,int(y_divisions/2)])
print ("Time index = ", time_index)
pl.show()
pl.clf()



params = maxwell.fit(abs(initial_conditions_velocity_x), floc=0)

pl.figure()
pl.hist(abs(initial_conditions_velocity_x), bins=30, normed=True)
initialConditionsVelocityX = np.linspace(0, 6, 100)
pl.plot(initialConditionsVelocityX, maxwell.pdf(initialConditionsVelocityX, *params), lw=3)
pl.show()

pl.clf()
 
