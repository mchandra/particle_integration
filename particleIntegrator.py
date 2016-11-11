import numpy as np
import pylab as pl
from scipy.integrate import odeint
import scipy.stats as stats
import matplotlib.pyplot as plt
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
""" Setting number of particles"""

no_of_particles = int(input("Enter number of particles:"))
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
x_divisions=int(input("Enter the number of divisions along x:"))
dx = (length_of_box_x/x_divisions)
x = np.arange(left_boundary,right_boundary,dx)
x= np.concatenate((x,[right_boundary]),axis = 0)

y_divisions=int(input("Enter the number of divisions along y:"))
dy = (length_of_box_y/y_divisions)
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
""" Initializing density matrix """

density = np.zeros((x.size-1,y.size-1),dtype=np.float)

#%%

"""Computing the density matrix for the initial conditions"""

i=0
j=0
x_zone = 0
y_zone = 0



for p in range(0,no_of_particles):
    for i in range(0,x.size-1):
        if((initial_conditions[p]>x[i])and(initial_conditions[p]<x[i+1])):
            x_zone = i
    for j in range(0,y.size-1):
        if((initial_conditions[p+no_of_particles]>y[j])and(initial_conditions[p+no_of_particles]<y[j+1])):
            y_zone = j
    
    density[x_zone,y_zone] = density[x_zone,y_zone] +1
       
        
       
density=density.transpose()
        
        
for i in range(0,y.size-1):
    temp_array=density[:,i]
    density[:,i] = temp_array[::-1]

#%%
""" normalizing the density """
density = density/(no_of_particles/(x_divisions*y_divisions))

""" sin perturbation """
for i in range(0,x.size-1):
    density[:,i] = density[:,i] + 0.5*np.sin((2*i*np.pi)/x_divisions)

"""computation of error along the centerline"""

error=np.zeros(density.size,dtype=np.float)
for i in range(0,density.size):
	error[i]=density[(y_divisions/2),i]-0.5*np.sin((2*i*np.pi)/x_divisions)-1



"""computation of standard deviation"""
total=0
for i in range(0,density.size):
	total=total+error[i]
mean=total/(density.size)
variance=0
for i in range(0,error.size):
	variance=variance+(error[i]-mean)**2
variance=variance/(error.size-1)
std_dev=math.sqrt(variance)
total=0
for i in range(0,density.size):
	total=total+abs(error[i])
mean_error=total/(density.size)

print("The variance of error is:",variance)
print("The standard deviation of error is:",std_dev)
print("The mean value of error is:",mean_error)


""" plotting the density contour plot"""
X, Y = np.meshgrid(x,y)


#pl.contourf(density,100)
#pl.colorbar()
#pl.savefig('/media/tejas/games/quazarWork/nParticleDensity/contour/point_mass' + '%04d'%time_index + '.png')
#pl.savefig('/media/tejas/games/quazarWork/nParticleDensity/InitialContour/IC'+'%04d'%x_divisions  + 'divisions.png')
#pl.show()

pl.figure()
pl.title('Numerical Density Along Centerline')
pl.xlabel('$x$')
pl.ylabel('$\mathrm{Density}$')
pl.plot(density[(y_divisions/2),:])
pl.show()

pl.figure()
pl.plot(error,label='Error')
pl.title('Error in Numerical Density Along Centerline')
pl.xlabel('$x$')
pl.ylabel('Error')
pl.show()

