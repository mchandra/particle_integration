import numpy as np
import h5py


""" Setting number of particles and other parameters"""

no_of_particles = 10000
x_coll = 100
y_coll = 1

x_divisions=100
y_divisions=1
z_divisions=1

""" Initial conditions """

left_boundary = 0
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x*np.random.rand(no_of_particles)

bottom_boundary = 0
top_boundary = 1
length_of_box_y           = length_of_box_x
initial_conditions_position_y = bottom_boundary + length_of_box_y*np.random.rand(no_of_particles)

back_boundary = 0
front_boundary = 1
length_of_box_z           = length_of_box_x
initial_conditions_position_z = back_boundary + length_of_box_z*np.random.rand(no_of_particles)


""" Setting velocities according to maxwellian distribution """

k=1
m=1
T=2
T_walls=2
const_gas=np.sqrt((k*T)/(m))
const_wall=np.sqrt((k*T_walls)/(m))
initial_conditions_velocity_x=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_y=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_z=np.zeros(no_of_particles,dtype=np.float)

for i in range(0,no_of_particles):
    r_1=np.random.rand(1)
    r_2=np.random.rand(1)
    r_3=np.random.rand(1)
    r_4=np.random.rand(1)
    
    y_1 = 4*r_1
    y_2 = 4*r_2
    y_3 = 4*r_3
    
    initial_conditions_velocity_x[i]=y_1
    initial_conditions_velocity_y[i]=y_2
    initial_conditions_velocity_z[i]=y_3




""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                                    initial_conditions_position_z, initial_conditions_velocity_x,\
                                    initial_conditions_velocity_y,initial_conditions_velocity_z],axis = 0)


h5f = h5py.File('initial_conditions.h5', 'w')
h5f.create_dataset('initial_conditions_dataset', data=initial_conditions)
h5f.close()


""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)
final_time            = 10 * box_crossing_time_scale
dt   = 0.01 * box_crossing_time_scale
time = np.arange(0, final_time, dt)  

""" Verlet Integrator """


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



x_div = 0
y_div = 0


""" Solving with collisions """

old= np.zeros(6*no_of_particles,dtype=np.float)



for time_index,t0 in enumerate(time):
    print("Computing for TimeIndex = ",time_index)
    t0 = time[time_index]
    list_indices = [[ [] for j in range(y_coll)  ] for i in range(x_coll)     ]
    if(time_index==time.size-1):
        break
    t1 = time[time_index+1]
    t = [t0, t1]
    if(time_index==0):
        initial_conditions = initial_conditions
    else:
        initial_conditions = old

    sol = Verlet(initial_conditions,t)


    for i in range(0,3*no_of_particles):
        if(sol[i]>=right_boundary):
            sol[i] = sol[i] - length_of_box_x
        if(sol[i]<=left_boundary):
            sol[i] = sol[i] + length_of_box_x
    

    
    for i in range(no_of_particles):
        x_div = int(100*sol[i])
        y_div = int(0.5)
            
        list_indices[x_div][y_div].append(i)

    for i in range(x_coll):
        for j in range(y_coll):
            if(len(list_indices[i][j])>1):
                for k in (list_indices[i][j]):
                    r_1=np.random.rand(1)
                    r_2=np.random.rand(1)
                    r_3=np.random.rand(1)
                    r_4=np.random.rand(1)
                    y_1=const_gas*np.sqrt(-2*np.log(r_1))*np.cos(2*np.pi*r_2)
                    y_2=const_gas*np.sqrt(-2*np.log(r_1))*np.sin(2*np.pi*r_2)
                    y_3=const_gas*np.sqrt(-2*np.log(r_3))*np.cos(2*np.pi*r_4)
                    sol[3*no_of_particles+k]=y_1
                    sol[4*no_of_particles+k]=y_2
                    sol[5*no_of_particles+k]=y_3


    old=sol
    list_indices.clear()
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=sol)
    h5f.close()



