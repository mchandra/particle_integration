import numpy as np
import pylab as pl
import scipy.stats as stats
import h5py
import math

""" Setting number of particles and other parameters"""

no_of_particles = int(input("Enter number of particles:"))
x_divisions=int(input("Enter the number of divisions along x:"))
y_divisions=int(input("Enter the number of divisions along y:"))

"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,t):
	x=initial_conditions[0:no_of_particles]
	y=initial_conditions[no_of_particles:2*no_of_particles]
	v_x=initial_conditions[2*no_of_particles:3*no_of_particles]
	v_y=initial_conditions[3*no_of_particles:4*no_of_particles]
	x_new = x + v_x*(t[1]-t[0])
	v_x_new = v_x # + -1*x_new*(t[1]-t[0]) #For Harmonic Oscillator
	y_new = y + v_y*(t[1]-t[0])
	v_y_new = v_y #+ -1*y_new*(t[1]-t[0]) #For Harmonic Oscillator
	nextstep=np.concatenate([x_new, y_new,v_x_new, v_y_new],axis=0)
	return(nextstep)



""" Initial conditions """

left_boundary = 0
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x=np.zeros(no_of_particles)
last=0
next=0
for i in range(x_divisions):
    next=last+(no_of_particles*0.5*np.sin(2*i*np.pi/x_divisions)/x_divisions)+(no_of_particles/x_divisions)
    initial_conditions_position_x[int(last):int(next)+1] = length_of_box_x*(i+1)/(x_divisions+1)
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




""" Discretizing time and making sure scaling is done right """
box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)
final_time            = 10 * box_crossing_time_scale
dt   = 0.01 * box_crossing_time_scale
time = np.arange(0, final_time, dt)



""" Initialization of solution matrix for all particles """


solution_all_current = np.zeros((4*no_of_particles),dtype=np.float)

#%%

""" Solving """
BC=input("1 for Hard and 2 for Periodic:")


for time_index,t0 in enumerate(time):
    
    t0 = time[time_index] 
    
    if(time_index==time.size-1):
        break
    t1 = time[time_index+1]
    t = [t0, t1] 
    if(time_index==0):
        initial_conditions = initial_conditions
       
    else:
        initial_conditions = (solution_all_current[:])
    sol = Verlet(initial_conditions,t)
    solution_all_current[:] = sol
    if(BC=='1'):
    
        for i in range(2*no_of_particles):
            if(solution_all_current[i]>right_boundary):
                solution_all_current[2*no_of_particles+i] = solution_all_current[2*no_of_particles+i] * (-1)
            if(solution_all_current[i,time_index]<left_boundary):
                solution_all_current[2*no_of_particles+i] = solution_all_current[2*no_of_particles+i] * (-1)
                
    if(BC=='2'):
        for i in range(2*no_of_particles):
            if(solution_all_current[i]>right_boundary):
                solution_all_current[i] = solution_all_current[i] - length_of_box_x
            if(solution_all_current[i]<left_boundary):
                solution_all_current[i] = solution_all_current[i] + length_of_box_x
                




    h5f = h5py.File('solution_all/solution'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=solution_all_current)
    h5f.close()
    
