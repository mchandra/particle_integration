import numpy as np
import scipy 
import random
import pylab as pl
from scipy.integrate import odeint
import scipy.integrate as integ
import scipy.stats as stats
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
import plotly.plotly as py

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
""" Setting number of particles"""

no_of_particles = 5000

#%%

""" function definition """

def dY_dt(Y, t):
    x = Y[0:no_of_particles]
    y = Y[no_of_particles:2*no_of_particles]
    dx_dt = Y[2*no_of_particles:3*no_of_particles]
    dy_dt = Y[3*no_of_particles:4*no_of_particles]
    dv_x_dt = -0*np.zeros(no_of_particles)
    dv_y_dt = -0*np.zeros(no_of_particles)
    return np.concatenate([dx_dt,dy_dt,dv_x_dt,dv_y_dt],axis =0)



#%%


""" Initial conditions """

left_boundary = -3
right_boundary = 3
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x*np.random.rand(no_of_particles)


bottom_boundary = -3
top_boundary = 3
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
""" Discretizing Space """

dx = 1
x = np.arange(left_boundary,right_boundary,dx)
x= np.concatenate((x,[right_boundary]),axis = 0)

dy = 1
y = np.arange(bottom_boundary,top_boundary,dy)
y = np.concatenate((y,[top_boundary]), axis =0)

#%%
""" Discretizing time and making sure scaling is done right """
box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)
final_time            = 5 * box_crossing_time_scale
dt   = 0.01 * box_crossing_time_scale
time = np.arange(0, final_time, dt)

a = np.zeros((4,6), dtype =np.float)


#%%

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
        print('started')
    else:
        initial_conditions = (solution_all[:,time_index-1])
    sol = odeint(dY_dt,initial_conditions,t)
    solution_all[:,time_index] = sol[1,:]
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
                


                

             
            
#%%
""" Initializing density matrix """

density = np.zeros((x.size-1,y.size-1,time.size),dtype=np.float)




#%%
""" Computing density plots """
i=0
j=0
x_zone = 0
y_zone = 0


for time_index,t0 in enumerate(time):
    for p in range(0,no_of_particles):
        for i in range(0,x.size-1):
            if((solution_all[p,time_index]>x[i])and(solution_all[p,time_index]<x[i+1])):
                x_zone = i
        for j in range(0,y.size-1):
            if((solution_all[p+no_of_particles,time_index]>y[j])and(solution_all[p+no_of_particles,time_index]<y[j+1])):
                y_zone = j
        
        density[x_zone,y_zone,time_index] = density[x_zone,y_zone,time_index] +1
       
        
       
for time_index,t0 in enumerate(time):
    density[:,:,time_index]= density[:,:,time_index].transpose()
        
        
for time_index,t0 in enumerate(time):
    for i in range(0,y.size-1):
        temp_array=density[:,i,time_index]
        density[:,i,time_index] = temp_array[::-1]


#%%

""" plotting the density contour plot"""
X, Y = np.meshgrid(x,y)

density_at_timeindex = density[:,:,0]
pl.contourf(density[:,:,0],100)
pl.colorbar()
pl.show()
#pl.savefig('/media/tejas/games/quazarWork/nParticleDensity/contour/point_mass' + '%04d'%time_index + '.png')
pl.clf()

#%%
""" Making images of the particles """


for time_index,t0 in enumerate(time[::10]):
    for i in range(0,no_of_particles):
        
        plt.plot(solution_all[i][time_index],solution_all[no_of_particles + i][time_index], 'o',color='blue', markersize=10, alpha = 0.4)
        plt.axhline(y=-1)
        plt.axvline(x=-1)
        plt.axhline(y=1)
        plt.axvline(x=1)
        plt.xlim([left_boundary, right_boundary])
        plt.ylim([bottom_boundary, top_boundary])
        plt.title('$\mathrm{Time}$ = ' + str(time[time_index]) )
        plt.xlabel('$x$')
        plt.ylabel('$y$')
    print ("Time index = ", time_index)
    plt.savefig('/media/tejas/games/quazarWork/nParticleDensity/images/point_mass' + '%04d'%time_index + '.png')
    plt.clf()
   
    
#%%

""" plotting maxwellian distribution curves """




params = maxwell.fit(abs(initial_conditions_velocity_x), floc=0)


plt.hist(abs(initial_conditions_velocity_x), bins=30, normed=True)
initialConditionsVelocityX = np.linspace(0, 6, 100)
plt.plot(initialConditionsVelocityX, maxwell.pdf(initialConditionsVelocityX, *params), lw=3)
plt.show()
#plt.savefig('/media/tejas/games/quazarWork/nParticleDensity/velocity_x_PDC')
plt.clf()

maxwell = stats.maxwell


params = maxwell.fit(abs(initial_conditions_velocity_y), floc=0)


plt.hist(abs(initial_conditions_velocity_y), bins=30, normed=True)
initialConditionsVelocityY = np.linspace(0, 6, 100)
plt.plot(initialConditionsVelocityY, maxwell.pdf(initialConditionsVelocityY, *params), lw=3)
plt.show()
#plt.savefig('/media/tejas/games/quazarWork/nParticleDensity/velocity_y_PDC')
plt.clf()





#%%
