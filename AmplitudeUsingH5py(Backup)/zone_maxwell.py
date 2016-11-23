from initialization import *
import numpy as np
import pylab as pl
import h5py
import scipy.stats as stats


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



""" Initializing """

particle_x_zone = np.zeros((no_of_particles),dtype=np.float)
particle_y_zone = np.zeros((no_of_particles),dtype=np.float)

""" Computing density plots """

i=0
j=0
x_zone = 0
y_zone = 0

""" Getting the x_zone and y_zone for each particle"""


for p in range(0,no_of_particles):
    print(" getting particle index"+str(p))
    for i in range(0,x.size-1):
        h5f = h5py.File('solution_all/solution_0.h5', 'r')
        solution_all_current = h5f['solution_all/solution_dataset_0'][:]
        h5f.close()
        if((solution_all_current[p]>x[i])and(solution_all_current[p]<x[i+1])):
            particle_x_zone[p] = i 
    for j in range(0,y.size-1):
        if((solution_all_current[p+no_of_particles]>y[j])and(solution_all_current[p+no_of_particles]<y[j+1])):
            particle_y_zone[p] = j


data_x_zone = []

for l in range(x.size-1):

    for p in range(0,no_of_particles):
        if(particle_x_zone[p]==l):
            data_x_zone.append(p)
    data_x_zone = np.array(data_x_zone)
    h5f = h5py.File('maxwell/max_zone_maxwell_'+str(l)+'.h5', 'w')
    h5f.create_dataset('maxwell/x_zone_maxwell_dataset_'+str(l), data=data_x_zone)
    h5f.close()
    data_x_zone.clear()

    
    
    
"""Plotting maxwell distribution for each zone """

maxwell_plot = []

for l in range(x.size-1):
    h5f = h5py.File('maxwell/max_zone_maxwell_'+str(l)+'.h5', 'r')
    data_x_zone = h5f['maxwell/x_zone_maxwell_dataset_'+str(l)][:]
    h5f.close()
    print("zone number : " + str(l))
    for index,i in range(data_x_zone):
        temp_velocity = np.zeros((data_x_zone.size),dtype=np.float)
        temp_velocity[index]=initial_conditions_velocity_y[i]
    
    
    
    
    params = maxwell.fit(abs(temp_velocity), floc=0)
    
    pl.figure()
    pl.hist(abs(temp_velocity), bins=30, normed=True)
    temp_velocity = np.linspace(0, 6, 100)
    pl.plot(temp_velocity, maxwell.pdf(temp_velocity, *params), lw=3)
    pl.savefig('maxwell/images/%04d'%l + '.png')

    pl.clf()


"""Plotting of the temperature distribution in each zome"""
