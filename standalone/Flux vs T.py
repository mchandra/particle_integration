import numpy as np
import h5py
from scipy.special import erfinv

""" Setting number of particles and other parameters"""

T_right = 1.1

no_of_particles = 100000
T_left          = 1
T_initial       = (T_right+1)*0.5 


""" Initial conditions """

left_boundary                 = 0
right_boundary                = 1
length_of_box_x               = right_boundary - left_boundary
initial_conditions_position_x = np.random.rand(no_of_particles)

bottom_boundary               = 0
top_boundary                  = 1
length_of_box_y               = top_boundary - bottom_boundary
initial_conditions_position_y = np.random.rand(no_of_particles)

""" Setting velocities according to maxwellian distribution """

x1=np.random.rand(no_of_particles)
x2=np.random.rand(no_of_particles) 

initial_conditions_velocity_x=np.sqrt(2*T_initial)*erfinv(2*x1-1)
initial_conditions_velocity_y=np.sqrt(2*T_initial)*erfinv(2*x2-1)

""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                                    initial_conditions_velocity_x, initial_conditions_velocity_y
                                    ],axis=0)

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x/np.max(initial_conditions_velocity_x)

final_time            = 10*box_crossing_time_scale
dt                    = 0.01*box_crossing_time_scale
time                  = np.arange(0, final_time, dt)

mom_x    = np.zeros(time.size)
mom_y    = np.zeros(time.size)
Energy   = np.zeros(time.size)
Pressure = np.zeros(time.size)
HeatFlux = np.zeros(time.size)


"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,dt):

  x = initial_conditions[0:no_of_particles]
  y = initial_conditions[no_of_particles:2*no_of_particles]
  
  v_x = initial_conditions[2*no_of_particles:3*no_of_particles]
  v_y = initial_conditions[3*no_of_particles:4*no_of_particles]
  
  x_new   = x + v_x*dt
  v_x_new = v_x
  y_new   = y + v_y*dt
  v_y_new = v_y
  
  nextstep = np.concatenate([x_new,y_new,v_x,v_y],axis=0)
  return(nextstep)

""" Initialization of solution matrix for all particles """

for time_index,t0 in enumerate(time):

  print("Computing for TimeIndex = ",time_index)

  if(time_index==time.size-1):
      break

  if(time_index==0):
      initial_conditions = initial_conditions
  else:
      initial_conditions = old

  sol = Verlet(initial_conditions,dt)

  x_coords = sol[0:no_of_particles]
  y_coords = sol[no_of_particles:2*no_of_particles]
  v_x      = sol[2*no_of_particles:3*no_of_particles]
  v_y      = sol[3*no_of_particles:4*no_of_particles]
  
  wall_x_left  = np.where(x_coords<left_boundary)[0]
  wall_x_right = np.where(x_coords>right_boundary)[0]
  wall_y_bot   = np.where(y_coords<bottom_boundary)[0]
  wall_y_top   = np.where(y_coords>top_boundary)[0]

  n1  = wall_x_left.size
  n2  = wall_x_right.size

  x1l = np.random.rand(n1) 
  x2l = np.random.rand(n1)
  x3l = np.random.rand(n1)
  
  x1r = np.random.rand(n2)
  x2r = np.random.rand(n2)
  x3r = np.random.rand(n2)


  sol[wall_x_left]  = 0
  sol[wall_x_right] = 1 - 1e-12
  
  sol[wall_x_left + 2*no_of_particles] = np.sqrt(-2*T_left*np.log(x1l))
  sol[wall_x_right+ 2*no_of_particles] = np.sqrt(-2*T_right*np.log(x1r))*(-1)

  sol[wall_x_left + 3*no_of_particles] = np.sqrt(2*T_left)*erfinv(2*x2l-1)
  sol[wall_x_right+ 3*no_of_particles] = np.sqrt(2*T_right)*erfinv(2*x2r-1)

  sol[wall_y_bot]   = sol[wall_y_bot] + 1
  sol[wall_y_top]   = sol[wall_y_top] - 1
    
  for i in range((no_of_particles-1)):
    j    = sol[(i+1):no_of_particles]        
    k    = sol[(i+1+no_of_particles):2*no_of_particles]
    velx = sol[i+2*no_of_particles]*np.ones(no_of_particles-1-i)
    vely = sol[i+3*no_of_particles]*np.ones(no_of_particles-1-i)
    
    velx_others = sol[(i+1+2*no_of_particles):3*no_of_particles]
    vely_others = sol[(i+1+3*no_of_particles):4*no_of_particles]
    
    x_particle = sol[i]*np.ones(j.size)
    y_particle = sol[i+no_of_particles]*np.ones(k.size)        
    j          = j-x_particle
    k          = k-y_particle
    dist       = np.sqrt(j**2+k**2)
    j          = j/dist
    k          = k/dist
    
    test_collision = dist<0.01 
    indices        = np.where(test_collision)[0]
    
    if(np.sum(test_collision)!=0):
      p           = (velx*j+vely*k-velx_others*j-vely_others*k)*test_collision
      velx        = velx*test_collision-p*j
      velx_others = velx_others+p*j
      vely        = vely*test_collision-p*k
      vely_others = vely_others+p*k
      index       = np.random.randint(0,indices.size)

      sol[i+1+indices[index]+2*no_of_particles] = velx_others[indices[index]]
      sol[i+1+indices[index]+3*no_of_particles] = vely_others[indices[index]]

      sol[i+2*no_of_particles] = velx[indices[index]]
      sol[i+3*no_of_particles] = vely[indices[index]]

  old = sol
  v_x = np.sum(sol[2*no_of_particles:3*no_of_particles])
  v_y = np.sum(sol[3*no_of_particles:4*no_of_particles])
  En  = np.sum(sol[2*no_of_particles:3*no_of_particles]**2 + sol[3*no_of_particles:4*no_of_particles]**2)
  qx  = np.sum(sol[2*no_of_particles:3*no_of_particles]*(sol[2*no_of_particles:3*no_of_particles]**2 + sol[3*no_of_particles:4*no_of_particles]**2))

  mom_x[time_index]    = v_x
  mom_y[time_index]    = v_y
  Energy[time_index]   = En/2
  Pressure[time_index] = En/no_of_particles
  HeatFlux[time_index] = qx/no_of_particles
  print("Pressure = ",Pressure[time_index])

  h5f = h5py.File('post/post_'+str(time_index)+'.h5', 'w')
  h5f.create_dataset('sol',     data=sol)
  h5f.create_dataset('mom_x',   data=mom_x)
  h5f.create_dataset('mom_y',   data=mom_y)
  h5f.create_dataset('Energy',  data=Energy)
  h5f.create_dataset('Pressure',data=Pressure)
  h5f.create_dataset('HeatFlux',data=HeatFlux)
  h5f.create_dataset('time',    data=time)
  h5f.close()