import numpy as np
import h5py


""" Setting number of particles and other parameters"""

no_of_particles = 1000000
x_divisions=32
y_divisions=1
z_divisions=1

""" Initial conditions """

left_boundary = 0
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x*np.random.rand(no_of_particles)
#last=0
#next=0

#for i in range(x_divisions):
#    next=last+(no_of_particles*0*np.sin(2*i*np.pi/x_divisions)/x_divisions)+(no_of_particles/x_divisions)
#    initial_conditions_position_x[int(round(last)):int(round(next))] = length_of_box_x*(2*i+1)/(2*x_divisions)
#    last=next

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
const=np.sqrt((k*T)/(m))
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

v_x=np.zeros(no_of_particles,dtype=np.float)
v_x=initial_conditions[3*no_of_particles:4*no_of_particles]

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(v_x)
final_time            = 20 * box_crossing_time_scale
dt   = 0.004 * box_crossing_time_scale
time = np.arange(0, final_time, dt)   

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

n = np.zeros(x.size-1,dtype=np.float)

""" Initialization of solution matrix for all particles """

sol = np.zeros(6*no_of_particles,dtype=np.float)

pressuredata = np.zeros(time.size,dtype=np.float)
heatfluxdata_x = np.zeros(time.size,dtype=np.float)
heatfluxdata_y = np.zeros(time.size,dtype=np.float)
heatfluxdata_z = np.zeros(time.size,dtype=np.float)


""" Solving """

old= np.zeros(6*no_of_particles,dtype=np.float)
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

    sol = Verlet(initial_conditions,t)
   
    
    for i in range(no_of_particles):
        alternator=0
        alternatol=0
        if(sol[i]>=right_boundary):#Cold Resevoir
            if(alternator%2==0):
                x_1=np.random.rand(1)
                x_2=np.random.rand(1)
                sol[3*no_of_particles+i] = abs(np.sqrt(1)*np.sqrt(-2*np.log(x_1))*np.cos(2*np.pi*x_2)) * (-1)
            else:
                x_1=np.random.rand(1)
                x_2=np.random.rand(1)
                sol[3*no_of_particles+i] = abs(np.sqrt(1)*np.sqrt(-2*np.log(x_1))*np.sin(2*np.pi*x_2)) * (-1)
            alternator=alternator+1
        if(sol[i]<=left_boundary):#Hot Resevoir
            if(alternatol%2==0):
                x_1=np.random.rand(1)
                x_2=np.random.rand(1)
                sol[3*no_of_particles+i] = abs(np.sqrt(3)*np.sqrt(-2*np.log(x_1))*np.cos(2*np.pi*x_2)) * (+1)
            else:
                x_1=np.random.rand(1)
                x_2=np.random.rand(1)
                sol[3*no_of_particles+i] = abs(np.sqrt(3)*np.sqrt(-2*np.log(x_1))*np.sin(2*np.pi*x_2)) * (+1)
            alternatol=alternatol+1
   
    for i in range(2*no_of_particles,3*no_of_particles):
        if(sol[i]>=right_boundary):
            sol[i] = sol[i] - length_of_box_x
        if(sol[i]<=left_boundary):
            sol[i] = sol[i] + length_of_box_x
    
    old=sol
    pressure=0
    heatflux_x=0
    heatflux_y=0
    heatflux_z=0
    for i in range(no_of_particles):
        pressure=pressure+sol[i+3*no_of_particles]**2+sol[i+4*no_of_particles]**2\
                                                     +sol[i+5*no_of_particles]**2
        heatflux_x=heatflux_x+sol[i+3*no_of_particles]*(sol[i+3*no_of_particles]**2+\
                 sol[i+4*no_of_particles]**2+sol[i+5*no_of_particles]**2)    
        heatflux_y=heatflux_y+sol[i+4*no_of_particles]*(sol[i+3*no_of_particles]**2+\
                 sol[i+4*no_of_particles]**2+sol[i+5*no_of_particles]**2)    
        heatflux_z=heatflux_z+sol[i+5*no_of_particles]*(sol[i+3*no_of_particles]**2+\
                 sol[i+4*no_of_particles]**2+sol[i+5*no_of_particles]**2)    
    heatflux=heatflux/no_of_particles
    pressure=pressure/no_of_particles
    print("Pressure = ",pressure)
    pressuredata[time_index]=pressure
    heatfluxdata_x[time_index]=heatflux_x
    heatfluxdata_y[time_index]=heatflux_y
    heatfluxdata_z[time_index]=heatflux_z
    print("Heat Flux in x-direction= ",heatflux_x)
        

h5f = h5py.File('post.h5', 'w')
h5f.create_dataset('time', data=time)
h5f.create_dataset('heatflux_x', data=heatfluxdata_x)
h5f.create_dataset('heatflux_y', data=heatfluxdata_y)
h5f.create_dataset('heatflux_z', data=heatfluxdata_z)
h5f.create_dataset('pressure', data=pressuredata)
h5f.close()
