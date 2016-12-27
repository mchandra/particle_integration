import numpy as np
import h5py


""" Setting number of particles and other parameters"""

no_of_particles = 10


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

const_gas=np.sqrt((k*T)/(m))

initial_conditions_velocity_x=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_y=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_z=np.zeros(no_of_particles,dtype=np.float)

for i in range(0,no_of_particles):
    
    
    
    r_1=np.random.rand(1)
    r_2=np.random.rand(1)
    r_3=np.random.rand(1)
    r_4=np.random.rand(1)
    
    y_1 = const_gas*np.sqrt(-2*np.log(r_1))*np.cos(2*np.pi*r_2)
    y_2 = const_gas*np.sqrt(-2*np.log(r_1))*np.sin(2*np.pi*r_2)
    y_3 = const_gas*np.sqrt(-2*np.log(r_3))*np.cos(2*np.pi*r_4)
    
    initial_conditions_velocity_x[i]=y_1
    initial_conditions_velocity_y[i]=y_2
    initial_conditions_velocity_z[i]=y_3
    #initial_conditions_velocity_x[i]=2
    #initial_conditions_velocity_y[i]=2
    #initial_conditions_velocity_z[i]=0




""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                                    initial_conditions_position_z, initial_conditions_velocity_x,\
                                    initial_conditions_velocity_y,initial_conditions_velocity_z],axis = 0)

""" Electric and Magnetic field """

Bx = 20
By = 0
Bz = 0
Ex = 0
Ey = 0
Ez = 0
me = 1
charge = 1
# mass of electron


""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)

final_time            = 30 * box_crossing_time_scale
dt   = 0.01 * box_crossing_time_scale
time = np.arange(0, final_time, dt)  

""" Boris method modified Verlet Integrator """


def mag_Verlet(initial_conditions,t):
    
    x=initial_conditions[0:no_of_particles]
    y=initial_conditions[no_of_particles:2*no_of_particles]
    z=initial_conditions[2*no_of_particles:3*no_of_particles]
    v_x=initial_conditions[3*no_of_particles:4*no_of_particles]
    v_y=initial_conditions[4*no_of_particles:5*no_of_particles]
    v_z=initial_conditions[5*no_of_particles:6*no_of_particles]
    x_new = x + v_x*(t[1]-t[0])
    y_new = y + v_y*(t[1]-t[0])
    z_new = z + v_z*(t[1]-t[0])    
    
    
    
    v_x_minus = v_x + (charge*Ex*(t[1]-t[0]))/(2*me)
    v_y_minus = v_y + (charge*Ey*(t[1]-t[0]))/(2*me)
    v_z_minus = v_z + (charge*Ez*(t[1]-t[0]))/(2*me)
    
    t_magx = (charge*Bx*(t[1]-t[0]))/(2*me)
    t_magy = (charge*By*(t[1]-t[0]))/(2*me)
    t_magz = (charge*Bz*(t[1]-t[0]))/(2*me)
   
   
    vminus_cross_t_x = (v_y_minus*t_magz)-(v_z_minus*t_magy)
    vminus_cross_t_y = -((v_x_minus*t_magz)-(v_z_minus*t_magx))
    vminus_cross_t_z = (v_x_minus*t_magy)-(v_y_minus*t_magx)
    
    v_dashx = v_x_minus + vminus_cross_t_x
    v_dashy = v_y_minus + vminus_cross_t_y
    v_dashz = v_z_minus + vminus_cross_t_z
    
    t_mag = np.sqrt(t_magx**2+t_magy**2+t_magz**2)
    
    
    
    s_x = (2* t_magx)/(1+abs(t_mag**2))
    s_y = (2* t_magy)/(1+abs(t_mag**2))
    s_z = (2* t_magz)/(1+abs(t_mag**2))
    
   

    v_x_plus = v_x_minus + ((v_dashy*s_z)-(v_dashz*s_y))
    v_y_plus = v_y_minus + -((v_dashx*s_z)-(v_dashz*s_x))
    v_z_plus = v_z_minus + ((v_dashx*s_y)-(v_dashy*s_x))
    
    
    
    v_x_new = v_x_plus + (charge*Ex*(t[1]-t[0]))/(2*me)
    v_y_new = v_y_plus + (charge*Ey*(t[1]-t[0]))/(2*me)
    v_z_new = v_z_plus + (charge*Ez*(t[1]-t[0]))/(2*me)    
    

    nextstep=np.concatenate([x_new, y_new, z_new ,v_x_new, v_y_new, v_z_new],axis=0)
    return(nextstep)






""" Solving """

old= np.zeros(6*no_of_particles,dtype=np.float)



""" Solver """



for time_index,t0 in enumerate(time):
    print("Computing for TimeIndex = ",time_index)
    t0 = time[time_index]
    if(time_index==time.size-1):
        break
    t1 = time[time_index+1]
    t = [t0, t1]
    if(time_index==0):
        initial_conditions = initial_conditions
    else:
        initial_conditions = old

    sol = mag_Verlet(initial_conditions,t)


    for i in range(3*no_of_particles):
        if(sol[i]>=right_boundary):
            sol[i] = sol[i] - 1
        if(sol[i]<=left_boundary):
            sol[i] =  sol[i] +1
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=sol)
    h5f.close()
    
    temp = 0
    for some_indice in range(no_of_particles):
        temp += (sol[3*no_of_particles+some_indice]**2+sol[4*no_of_particles+some_indice]**2+sol[5*no_of_particles+some_indice]**2)/(3*no_of_particles)
        

    print('temperature at current timestep ', temp)



    old=sol


print(' the time scale is ', dt)
