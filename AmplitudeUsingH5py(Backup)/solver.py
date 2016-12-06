import numpy as np
import pylab as pl
import scipy.stats as stats
import h5py
from scipy.integrate import odeint

no_of_particles = 10000
x_divisions=32
y_divisions=1
length_of_box_x=1

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

"""Function for OdeInt"""

def dY_dt(Y, t):
    x = Y[0:no_of_particles]
    y = Y[no_of_particles:2*no_of_particles]
    dx_dt = Y[2*no_of_particles:3*no_of_particles]
    dy_dt = Y[3*no_of_particles:4*no_of_particles]
    dv_x_dt = -0*np.zeros(no_of_particles)
    dv_y_dt = -0*np.zeros(no_of_particles)
    return np.concatenate([dx_dt,dy_dt,dv_x_dt,dv_y_dt],axis =0)

""" Initialization of solution matrix for all particles """


solution_all_current = np.zeros((4*no_of_particles),dtype=np.float)

#%%


""" Solving """
BC='2'


h5f = h5py.File('initial_conditions.h5', 'r')
initial_conditions = h5f['initial_conditions_dataset'][:]
h5f.close()
print(initial_conditions,"	",initial_conditions.size)    
v_x=np.zeros(no_of_particles,dtype=np.float)
v_x=initial_conditions[2*no_of_particles:3*no_of_particles]

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(v_x)
final_time            = 5 * box_crossing_time_scale
dt   = 0.01 * box_crossing_time_scale
time = np.arange(0, final_time, dt)   
    
choice=input("Enter your choice: 1-Verlet, 2-Odeint:")

ic=np.zeros(4*no_of_particles,dtype=np.float)

for time_index,t0 in enumerate(time):
    
    t0 = time[time_index] 
    print("Computing for TimeIndex = ",time_index)
    if(time_index==time.size-1):
        break
    t1 = time[time_index+1]
    t = [t0, t1] 
    
    if(time_index==0):
        ic=initial_conditions
    else:
        ic=solution_all_current
       
    if(choice=='1'):
    
        sol = Verlet(ic,t)
        
        solution_all_current = sol
        if(BC=='1'):
        
            for i in range(2*no_of_particles):
                if(solution_all_current[i]>1):
                    solution_all_current[2*no_of_particles+i] = solution_all_current[2*no_of_particles+i] * (-1)
                if(solution_all_current[i]<0):
                    solution_all_current[2*no_of_particles+i] = solution_all_current[2*no_of_particles+i] * (-1)
                    
        if(BC=='2'):
            for i in range(2*no_of_particles):
                if(solution_all_current[i]>=1):
                    solution_all_current[i] = solution_all_current[i] - length_of_box_x
                if(solution_all_current[i]<0):
                    solution_all_current[i] = solution_all_current[i] + length_of_box_x
        
        
        
    if(choice=='2'):    

        sol = odeint(dY_dt,ic,t)
        solution_all_current = sol[1,:]
        if(BC=='1'):

            for i in range(2*no_of_particles):
                if(solution_all_current[i]>=1):
                    solution_all_current[2*no_of_particles+i] = solution_all_current[2*no_of_particles+i] * (-1)
                if(sol[i]<=0):
                    solution_all_current[2*no_of_particles+i] = solution_all_current[2*no_of_particles+i] * (-1)

        if(BC=='2'):
            for i in range(2*no_of_particles):
                if(solution_all_current[i]>=1):
                    print ("particle number = ", i, "  position right = ", solution_all_current[i])
                    solution_all_current[i] = solution_all_current[i] - length_of_box_x
                if(solution_all_current[i]<=0):
                    print ("particle number = ", i, "  position left = ", solution_all_current[i])
                    solution_all_current[i] = solution_all_current[i] + length_of_box_x

        



    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=solution_all_current)
    h5f.close()

