import numpy as np
import pylab as pl
import scipy.stats as stats
import math

#%%
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

#%%
""" Setting number of particles"""

no_of_particles = int(input("Enter number of particles:"))

#%%

"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,t):
	x=initial_conditions[0:no_of_particles]
	y=initial_conditions[no_of_particles:2*no_of_particles]
	v_x=initial_conditions[2*no_of_particles:3*no_of_particles]
	v_y=initial_conditions[3*no_of_particles:4*no_of_particles]
	x_new = x + v_x*(t[1]-t[0])
	v_x_new = v_x #+ -1*x_new*(t[1]-t[0]) For Harmonic Oscillator
	y_new = y + v_y*(t[1]-t[0])
	v_y_new = v_y #+ -1*y_new*(t[1]-t[0]) For Harmonic Oscillator
	nextstep=np.concatenate([x_new, y_new,v_x_new, v_y_new],axis=0)
	return(nextstep)
	
#%%


""" Initial conditions """

left_boundary = -1
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x*np.random.rand(no_of_particles)


bottom_boundary = -1
top_boundary = 1
length_of_box_y           = length_of_box_x
initial_conditions_position_y = bottom_boundary + length_of_box_y*np.random.rand(no_of_particles)


#%%

""" Setting velocities according to maxwellian distribution """
maxwell = stats.maxwell
maxwell.mean(loc=0, scale=1)

#%%
""" Setting the  velocity distribution"""

a=np.random.choice([1,-1], size = no_of_particles)
b=np.random.choice([1,-1], size = no_of_particles)
initial_conditions_velocity_x =  np.multiply(a,maxwell.rvs(loc=0, scale=1, size=no_of_particles))
initial_conditions_velocity_y =  np.multiply(b,maxwell.rvs(loc=0, scale=1, size=no_of_particles))


#%%
""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,\
                                    initial_conditions_position_y, initial_conditions_velocity_x,\
                                    initial_conditions_velocity_y],axis = 0)



#%%
""" Discretizing time and making sure scaling is done right """
box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)
final_time            = 10 * box_crossing_time_scale
dt   = 0.01 * box_crossing_time_scale
time = np.arange(0, final_time, dt)

""" Initialization of solution matrix for all particles """


solution_all = np.zeros((4*no_of_particles,time.size),dtype=np.float)

#%%

""" Solving """
BC=input("1 for Hard and 2 for Periodic")

for time_index,t0 in enumerate(time):
    
    t0 = time[time_index] 
    
    if(time_index==time.size-1):
        break
    t1 = time[time_index+1]
    t = [t0, t1] 
    if(time_index==0):
        initial_conditions = initial_conditions
       
    else:
        initial_conditions = (solution_all[:,time_index-1])
    sol = Verlet(initial_conditions,t)
    solution_all[:,time_index] = sol
    if(BC=='1'):
    
        for i in range(2*no_of_particles):
            if(solution_all[i,time_index]>right_boundary):
                solution_all[2*no_of_particles+i,time_index] = solution_all[2*no_of_particles+i,time_index] * (-1)
            if(solution_all[i,time_index]<left_boundary):
                solution_all[2*no_of_particles+i,time_index] = solution_all[2*no_of_particles+i,time_index] * (-1)
                
    if(BC=='2'):
        for i in range(2*no_of_particles):
            if(solution_all[i,time_index]>right_boundary):
                solution_all[i,time_index] = solution_all[i,time_index] - length_of_box_x
            if(solution_all[i,time_index]<left_boundary):
                solution_all[i,time_index] = solution_all[i,time_index] + length_of_box_x

"""Computation of Energy and Momentum"""
energy=np.zeros(time.size)
for time_index in enumerate(time):
    for i in range(no_of_particles):
        energy[time_index]=energy[time_index]+(solution_all[(i+2*no_of_particles),time_index])**2+(solution_all[(i+3*no_of_particles),time_index])**2
                      # +0*solution_all[i+no_of_particles,time_index]**2+0*solution_all[i,time_index]**2

mom=np.zeros(time.size)
for time_index in enumerate(time):
    for i in range(no_of_particles):
        mom[time_index]=mom[time_index]+math.sqrt((solution_all[(i+2*no_of_particles),time_index])**2+(solution_all[(i+3*no_of_particles),time_index])**2)

"""Plotting of the Energy and Momentum"""

pl.figure()
pl.plot(time,energy)
pl.title('Energy of Particle')
pl.xlabel('$\mathrm{Time}$')
pl.ylabel('$E$')
pl.show()

pl.figure()
pl.plot(time,mom)
pl.title('Momentum of Particle')
pl.xlabel('$\mathrm{Time}$')
pl.ylabel('$p$')
pl.show()

""" Making images of the particles """


for time_index,t0 in enumerate(time[::10]):
    for i in range(0,no_of_particles):
        
        pl.plot(solution_all[i][time_index],solution_all[no_of_particles + i][time_index], 'o',color='blue', markersize=10, alpha = 0.4)
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
    pl.savefig('point_mass' + '%04d'%time_index + '.png')
    pl.clf()






