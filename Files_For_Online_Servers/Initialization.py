import numpy as np
import h5py

""" Setting number of particles and other parameters"""

no_of_particles = 10000000
x_divisions=32
y_divisions=1
z_divisions=1

""" Initial conditions """

left_boundary = 0
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x=np.zeros(no_of_particles)
last=0
next=0

for i in range(x_divisions):
    next=last+(no_of_particles*0.5*np.sin(2*i*np.pi/x_divisions)/x_divisions)+(no_of_particles/x_divisions)
    initial_conditions_position_x[int(round(last)):int(round(next))] = length_of_box_x*(2*i+1)/(2*x_divisions)
    last=next

bottom_boundary = 0
top_boundary = 1
length_of_box_y           = length_of_box_x
initial_conditions_position_y = bottom_boundary + length_of_box_y*np.random.rand(no_of_particles)

back_boundary = 0
front_boundary = 1
length_of_box_z           = length_of_box_x
initial_conditions_position_z = back_boundary + length_of_box_z*np.random.rand(no_of_particles)

""" Discretizing Space """

dx = (length_of_box_x/x_divisions)
x = np.arange(left_boundary,right_boundary,dx)
x= np.concatenate((x,[right_boundary]),axis = 0)


""" Setting velocities according to maxwellian distribution """

k=1
m=1
T=2
const=np.sqrt((k*T)/(2*m))
initial_conditions_velocity_x=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_y=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_z=np.zeros(no_of_particles,dtype=np.float)

for i in range(0,no_of_particles,2):
    x_1=np.random.rand(1)
    x_2=np.random.rand(1)
    y_1=const*np.sqrt(-2*np.log(x_1))*np.cos(2*np.pi*x_2)
    y_2=const*np.sqrt(-2*np.log(x_1))*np.sin(2*np.pi*x_2)
    x_3=np.random.rand(1)
    x_4=np.random.rand(1)
    y_3=const*np.sqrt(-2*np.log(x_3))*np.cos(2*np.pi*x_4)
    y_4=const*np.sqrt(-2*np.log(x_3))*np.sin(2*np.pi*x_4)
    x_5=np.random.rand(1)
    x_6=np.random.rand(1)
    y_5=const*np.sqrt(-2*np.log(x_5))*np.cos(2*np.pi*x_6)
    y_6=const*np.sqrt(-2*np.log(x_5))*np.sin(2*np.pi*x_6)
    initial_conditions_velocity_x[i]=y_1
    initial_conditions_velocity_y[i]=y_2
    initial_conditions_velocity_z[i]=y_3
    initial_conditions_velocity_x[i+1]=y_4
    initial_conditions_velocity_y[i+1]=y_5
    initial_conditions_velocity_z[i+1]=y_6

""" Setting the  velocity distyribution"""

a=np.random.choice([1,-1], size = no_of_particles)
b=np.random.choice([1,-1], size = no_of_particles)
c=np.random.choice([1,-1], size = no_of_particles)
initial_conditions_velocity_x = np.multiply(initial_conditions_velocity_x,a) 
initial_conditions_velocity_y = np.multiply(initial_conditions_velocity_y,b) 
initial_conditions_velocity_z = np.multiply(initial_conditions_velocity_z,c) 
 

""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                                    initial_conditions_position_z, initial_conditions_velocity_x,\
                                    initial_conditions_velocity_y,initial_conditions_velocity_z],axis = 0)

""" printing out the initial conditions to a file """

h5f = h5py.File('initial_conditions.h5', 'w')
h5f.create_dataset('initial_conditions_dataset', data=initial_conditions)
h5f.close()
