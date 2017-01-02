import numpy as np
import numexpr as ne
import h5py
from scipy.special import erfinv
ne.set_num_threads(4)

""" Setting number of particles and other parameters"""
no_of_particles = 500000

""" Initial conditions """

left_boundary = 0
right_boundary = 1.0
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x*np.random.rand(no_of_particles)

bottom_boundary = 0
top_boundary = 1.0
length_of_box_y           = top_boundary - bottom_boundary
initial_conditions_position_y = bottom_boundary + length_of_box_y*np.random.rand(no_of_particles)

back_boundary = 0
front_boundary = 1.0
length_of_box_z           = front_boundary - back_boundary
initial_conditions_position_z = back_boundary + length_of_box_z*np.random.rand(no_of_particles)

""" Setting velocities according to maxwellian distribution """

k=1.0
m=1.0
T=1.5
T_walls=2.0

initial_conditions_velocity_x=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_y=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_z=np.zeros(no_of_particles,dtype=np.float)

x1 = np.random.rand(no_of_particles)
x2 = np.random.rand(no_of_particles)
x3 = np.random.rand(no_of_particles)
x4 = np.random.rand(no_of_particles)


initial_conditions_velocity_x=ne.evaluate("sqrt(2*T)*sqrt(-1*log(x2))*cos(2*3.14159*x1)")*np.random.choice([1,-1],size=no_of_particles)
initial_conditions_velocity_y=ne.evaluate("sqrt(2*T)*sqrt(-1*log(x2))*sin(2*3.14159*x1)")*np.random.choice([1,-1],size=no_of_particles)
initial_conditions_velocity_z=ne.evaluate("sqrt(2*T)*sqrt(-1*log(x4))*cos(2*3.14159*x3)")*np.random.choice([1,-1],size=no_of_particles)

""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                                    initial_conditions_position_z, initial_conditions_velocity_x,\
                                    initial_conditions_velocity_y,initial_conditions_velocity_z],axis = 0)


""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)
final_time            = 25
dt   = 0.05 * box_crossing_time_scale
time = np.arange(0, final_time, dt)
print(time.size)

"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,t):

    x=initial_conditions[0:no_of_particles]
    y=initial_conditions[no_of_particles:2*no_of_particles]
    z=initial_conditions[2*no_of_particles:3*no_of_particles]
    v_x=initial_conditions[3*no_of_particles:4*no_of_particles]
    v_y=initial_conditions[4*no_of_particles:5*no_of_particles]
    v_z=initial_conditions[5*no_of_particles:6*no_of_particles]
    dt=t[1]-t[0]
    x_new = ne.evaluate("x + v_x*dt")
    v_x_new = ne.evaluate("v_x")
    y_new = ne.evaluate("y + v_y*dt")
    v_y_new = ne.evaluate("v_y")
    z_new = ne.evaluate("z + v_z*dt")
    v_z_new = ne.evaluate("v_z")
    nextstep=np.concatenate([x_new, y_new, z_new ,v_x_new, v_y_new, v_z_new],axis=0)
    return(nextstep)

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

    x_coords=sol[0:no_of_particles]
    y_coords=sol[no_of_particles:2*no_of_particles]
    z_coords=sol[2*no_of_particles:3*no_of_particles]
    v_x=sol[3*no_of_particles:4*no_of_particles]
    v_y=sol[4*no_of_particles:5*no_of_particles]
    v_z=sol[5*no_of_particles:6*no_of_particles]

    """Implementing the Thermal B.C's"""
    test_left=ne.evaluate("x_coords<left_boundary")
    test_right=ne.evaluate("x_coords>right_boundary")

    test_left_complement=ne.evaluate("where(test_left,0,1)")

    test_right_complement=ne.evaluate("where(test_right,0,1)")
    test_left=ne.evaluate("where(test_left,1,0)")
    test_right=ne.evaluate("where(test_right,1,0)")

    collided=ne.evaluate("test_left_complement==test_right_complement")
    collided_complement=ne.evaluate("where(collided,0,1)")
    collided=ne.evaluate("where(collided,1,0)")
    
    x1 = np.random.rand(no_of_particles)
    x2 = np.random.rand(no_of_particles)
    x3 = np.random.rand(no_of_particles)
    x4 = np.random.rand(no_of_particles)

    #v_x=ne.evaluate("v_x*collided+test_left*sqrt(-2*T_walls*log(x1))-test_right*sqrt(-2*T_walls*log(x1))")

    v_x=ne.evaluate("v_x*collided+test_left*sqrt(2*T_walls)*sqrt(-1*log(x4))*cos(2*3.14159*x1)-\
                    test_right*sqrt(2*T_walls)*sqrt(-1*log(x4))*cos(2*3.14159*x1)")

    v_y=ne.evaluate("v_y*collided+test_left*sqrt(2*T_walls)*sqrt(-1*log(x2))*cos(2*3.14159*x3)+\
                   test_right*sqrt(2*T_walls)*sqrt(-1*log(x2))*cos(2*3.14159*x3)")*\
                    (collided_complement*np.random.choice([1,-1],size=no_of_particles)+collided)

    v_z=ne.evaluate("v_z*collided+test_left*sqrt(2*T_walls)*sqrt(-1*log(x2))*sin(2*3.14159*x3)+\
                    test_right*sqrt(2*T_walls)*sqrt(-1*log(x2))*sin(2*3.14159*x3)")*\
                   (collided_complement*np.random.choice([1,-1],size=no_of_particles)+collided)       

    sol[3*no_of_particles:4*no_of_particles]=v_x
    sol[4*no_of_particles:5*no_of_particles]=v_y
    sol[5*no_of_particles:6*no_of_particles]=v_z

    """Implementation of Periodic B.C's"""
  
    wall_y_bot=ne.evaluate("y_coords<bottom_boundary")
    wall_y_top=ne.evaluate("y_coords>top_boundary")
    wall_z_back=ne.evaluate("z_coords<back_boundary")
    wall_z_front=ne.evaluate("z_coords>front_boundary")

    wall_y_bot=ne.evaluate("where(wall_y_bot,1,0)")
    wall_y_top=ne.evaluate("where(wall_y_top,1,0)")
    wall_z_back=ne.evaluate("where(wall_z_back,1,0)")
    wall_z_front=ne.evaluate("where(wall_z_front,1,0)")

    y_coords=ne.evaluate("y_coords+wall_y_bot*length_of_box_y-wall_y_top*length_of_box_y")
    z_coords=ne.evaluate("z_coords+wall_z_back*length_of_box_z-wall_z_front*length_of_box_z")
    
    sol[no_of_particles:2*no_of_particles]=y_coords
    sol[2*no_of_particles:3*no_of_particles]=z_coords

    old=sol
 

    """Calculation of Pressure and Heat Fluxes"""

    pressure=0
    heatflux_x=0
    heatflux_y=0
    heatflux_z=0
             
    pressure=pressure+ne.evaluate("sum(v_x**2+v_y**2+v_z**2)")
    heatflux_x=heatflux_x+ne.evaluate("sum(v_x*(v_x**2+v_y**2+v_z**2))")
    heatflux_y=heatflux_y+ne.evaluate("sum(v_y*(v_x**2+v_y**2+v_z**2))") 
    heatflux_z=heatflux_z+ne.evaluate("sum(v_z*(v_x**2+v_y**2+v_z**2))")
        
    heatflux_x=heatflux_x/no_of_particles
    heatflux_y=heatflux_y/no_of_particles
    heatflux_z=heatflux_z/no_of_particles
    pressure=pressure/no_of_particles

    print("Pressure = ",pressure)
    
    pressuredata[time_index]=pressure
    heatfluxdata_x[time_index]=heatflux_x
    heatfluxdata_y[time_index]=heatflux_y
    heatfluxdata_z[time_index]=heatflux_z

    
h5f = h5py.File('post.h5', 'w')
h5f.create_dataset('time', data=time)
h5f.create_dataset('heatflux_x', data=heatfluxdata_x)
h5f.create_dataset('heatflux_y', data=heatfluxdata_y)
h5f.create_dataset('heatflux_z', data=heatfluxdata_z)
h5f.create_dataset('pressure', data=pressuredata)
h5f.close()
