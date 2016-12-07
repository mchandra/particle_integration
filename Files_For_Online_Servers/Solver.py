import numpy as np
import h5py

no_of_particles =1000000
x_divisions=32
y_divisions=1
length_of_box_x=1

"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,t):
    x=initial_conditions[0:no_of_particles]
    y=initial_conditions[no_of_particles:2*no_of_particles]
    z=initial_conditions[2*no_of_particles:3*no_of_particles]
    v_x=initial_conditions[3*no_of_particles:4*no_of_particles]
    v_y=initial_conditions[4*no_of_particles:5*no_of_particles]
    v_z=initial_conditions[5*no_of_particles:6*no_of_particles]
    x_new = x + v_x*(t[1]-t[0])
    v_x_new = v_x 
    y_new = y + v_y*(t[1]-t[0])
    v_y_new = v_y 
    z_new = z + v_z*(t[1]-t[0])
    v_z_new = v_z 
    nextstep=np.concatenate([x_new, y_new, z_new ,v_x_new, v_y_new, v_z_new],axis=0)
    return(nextstep)

""" Initialization of solution matrix for all particles """


solution_all_current = np.zeros((6*no_of_particles),dtype=np.float)
sol = np.zeros((6*no_of_particles),dtype=np.float)

""" Solving """
BC='2'


h5f = h5py.File('initial_conditions.h5', 'r')
initial_conditions = h5f['initial_conditions_dataset'][:]
h5f.close()
v_x=np.zeros(no_of_particles,dtype=np.float)
v_x=initial_conditions[3*no_of_particles:4*no_of_particles]

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(v_x)
final_time            = 4 * box_crossing_time_scale
dt   = 0.005 * box_crossing_time_scale
time = np.arange(0, final_time, dt)   
    
ic=np.zeros(6*no_of_particles,dtype=np.float)

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
       
         
    sol = Verlet(ic,t)
        
    solution_all_current = sol
    if(BC=='1'):
        
        for i in range(3*no_of_particles):
            if(solution_all_current[i]>1):
                solution_all_current[3*no_of_particles+i] = solution_all_current[3*no_of_particles+i] * (-1)
            if(solution_all_current[i]<0):
                solution_all_current[3*no_of_particles+i] = solution_all_current[3*no_of_particles+i] * (-1)
                    
    if(BC=='2'):
        for i in range(3*no_of_particles):
            if(solution_all_current[i]>=1):
                solution_all_current[i] = solution_all_current[i] - length_of_box_x
            if(solution_all_current[i]<0):
                solution_all_current[i] = solution_all_current[i] + length_of_box_x
    
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=solution_all_current)
    h5f.close()

