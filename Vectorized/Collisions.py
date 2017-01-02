import numpy as np
import h5py
from scipy.special import erfinv
import numexpr as ne
#import pylab as pl
ne.set_num_threads(4)

""" Setting number of particles and other parameters"""

no_of_particles = 1000

""" Initial conditions """

left_boundary = 0
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x = length_of_box_x*np.random.rand(no_of_particles)

bottom_boundary = 0
top_boundary = 1
length_of_box_y           = length_of_box_x
initial_conditions_position_y = length_of_box_y*np.random.rand(no_of_particles)

""" Setting velocities according to maxwellian distribution """

initial_conditions_velocity_x=np.random.rand(no_of_particles)*np.random.choice([-1,1],size=no_of_particles)
initial_conditions_velocity_y=np.random.rand(no_of_particles)*np.random.choice([-1,1],size=no_of_particles)

""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                                    initial_conditions_velocity_x,initial_conditions_velocity_y],axis = 0)

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x/np.max(initial_conditions_velocity_x)
final_time            = 20*box_crossing_time_scale
dt   = 0.0001*box_crossing_time_scale
time = np.arange(0, final_time, dt)

"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,t):
    x=initial_conditions[0:no_of_particles]
    y=initial_conditions[no_of_particles:2*no_of_particles]
    v_x=initial_conditions[2*no_of_particles:3*no_of_particles]
    v_y=initial_conditions[3*no_of_particles:4*no_of_particles]
    dt=t[1]-t[0]
    x_new = ne.evaluate("x + v_x*dt")
    v_x_new = ne.evaluate("v_x")
    y_new = ne.evaluate("y + v_y*dt")
    v_y_new = ne.evaluate("v_y")
    nextstep=np.concatenate([x_new, y_new ,v_x_new, v_y_new],axis=0)
    return(nextstep)

""" Initialization of solution matrix for all particles """

sol = np.zeros(4*no_of_particles,dtype=np.float)
fsol = np.zeros(4*no_of_particles,dtype=np.float)
""" Solving """

old= np.zeros(4*no_of_particles,dtype=np.float)
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
    v_x=sol[2*no_of_particles:3*no_of_particles]
    v_y=sol[3*no_of_particles:4*no_of_particles]
    
    wall_x_left=ne.evaluate("x_coords<left_boundary")
    wall_x_right=ne.evaluate("x_coords>right_boundary")
    wall_y_bot=ne.evaluate("y_coords<bottom_boundary")
    wall_y_top=ne.evaluate("y_coords>top_boundary")

    wall_x_left=ne.evaluate("where(wall_x_left,1,0)")
    wall_x_right=ne.evaluate("where(wall_x_right,1,0)")
    wall_y_bot=ne.evaluate("where(wall_y_bot,1,0)")
    wall_y_top=ne.evaluate("where(wall_y_top,1,0)")

    x_coords=ne.evaluate("x_coords+wall_x_left*length_of_box_x-wall_x_right*length_of_box_x")
    y_coords=ne.evaluate("y_coords+wall_y_bot*length_of_box_y-wall_y_top*length_of_box_y")
    
    
    sol[0:no_of_particles]=x_coords
    sol[no_of_particles:2*no_of_particles]=y_coords

    
    for i in range(no_of_particles):
        j=sol[(i+1):no_of_particles]        
        k=sol[(i+1+no_of_particles):2*no_of_particles]
        velx=sol[i+2*no_of_particles]*np.ones(j.size)
        vely=sol[i+3*no_of_particles]*np.ones(j.size)
        velx_others=sol[(i+1+2*no_of_particles):3*no_of_particles]
        vely_others=sol[(i+1+3*no_of_particles):4*no_of_particles]
        x_particle=sol[i]
        y_particle=sol[i+no_of_particles]        
        j=ne.evaluate("j-x_particle")
        k=ne.evaluate("k-y_particle")
        dist=ne.evaluate("sqrt(j**2+k**2)")
        j=j/dist
        k=k/dist
        test_collision=ne.evaluate("dist<0.02") 
        test_collision=ne.evaluate("where(test_collision,1,0)")
        p=(velx*j+vely*k-velx_others*j-vely_others*k)*test_collision
        velx=velx*test_collision-p*j
        velx_others=velx_others+p*j
        vely=vely*test_collision-p*k
        vely_others=vely_others+p*k 
        sol[(i+1+2*no_of_particles):3*no_of_particles]=velx_others
        sol[(i+1+3*no_of_particles):4*no_of_particles]=vely_others
        if(ne.evaluate("sum(test_collision)")==0):
            break
        else:
            sol[i+2*no_of_particles]=np.sum(velx)/np.sum(test_collision)
            sol[i+3*no_of_particles]=np.sum(vely)/np.sum(test_collision)

    old=sol
    """
    for i in range(no_of_particles):
        pl.plot(sol[i],sol[no_of_particles + i], 'o',color='blue', markersize=10, alpha = 0.4)
        pl.xlim([left_boundary, right_boundary])
        pl.ylim([bottom_boundary, top_boundary])
        pl.title('$\mathrm{Time}$ = ' + str(time[time_index]) )
        pl.xlabel('$x$')
        pl.ylabel('$y$')
    print ("Time index = ", time_index)
    pl.savefig('images/' + '%04d'%time_index + '.png')
    pl.clf() """

    if(time_index==time.size-2):
        fsol=sol    
h5f = h5py.File('post-collision.h5', 'w')
h5f.create_dataset('sol', data=fsol)
h5f.close()
