import numpy as np
import pylab as pl
import scipy.stats as stats
import h5py

""" Setting number of particles and other parameters"""

no_of_particles = 5000
x_divisions=32
y_divisions=1


""" Initial conditions """

left_boundary = 0
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x=np.zeros(no_of_particles)
last=0
next=0
for i in range(x_divisions):
    next=last+(no_of_particles*0.5*np.sin(2*i*np.pi/x_divisions)/x_divisions)+(no_of_particles/x_divisions)
    initial_conditions_position_x[int(last):int(next)] = length_of_box_x*(2*i+1)/(2*x_divisions)
    last=next

bottom_boundary = 0
top_boundary = 1
length_of_box_y           = length_of_box_x
initial_conditions_position_y = bottom_boundary + length_of_box_y*np.random.rand(no_of_particles)

""" Discretizing Space """

dx = (length_of_box_x/x_divisions)
x = np.arange(left_boundary,right_boundary,dx)
x= np.concatenate((x,[right_boundary]),axis = 0)
  

dy = (length_of_box_y/y_divisions)
y = np.arange(bottom_boundary,top_boundary,dy)
y = np.concatenate((y,[top_boundary]), axis =0)


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


""" printing out the initial conditions to a file """

h5f = h5py.File('initial_conditions.h5', 'w')
h5f.create_dataset('initial_conditions_dataset', data=initial_conditions)
h5f.close()
