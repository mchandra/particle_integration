import numpy as np
import h5py
import arrayfire as af


""" Setting number of particles and other parameters"""

no_of_particles = 10000

T_initial       = 1.5
T_left          = 1
T_right         = 2

""" Initial conditions """

left_boundary                 = 0
right_boundary                = 1
length_of_box_x               = right_boundary - left_boundary
initial_conditions_position_x = af.randu(no_of_particles)

bottom_boundary               = 0
top_boundary                  = 1
length_of_box_y               = top_boundary - bottom_boundary
initial_conditions_position_y = af.randu(no_of_particles)

""" Setting velocities according to maxwellian distribution """

x1=af.randu(no_of_particles)
x2=af.randu(no_of_particles) 

initial_conditions_velocity_x=np.sqrt(T_initial)*af.arith.sqrt(-2*af.arith.log(x1))*af.arith.cos(2*np.pi*x2)
initial_conditions_velocity_y=np.sqrt(T_initial)*af.arith.sqrt(-2*af.arith.log(x1))*af.arith.sin(2*np.pi*x2)

""" Combining the initial conditions into one vector"""

initial_conditions = af.join(0,initial_conditions_position_x,initial_conditions_position_y,\
                              initial_conditions_velocity_x, initial_conditions_velocity_y
                            )

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x/af.algorithm.max(initial_conditions_velocity_x)

final_time            = 10*box_crossing_time_scale
dt                    = 0.005*box_crossing_time_scale
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
  
  nextstep = af.join(0,x_new,y_new,v_x,v_y)
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
  
  wall_x_left  = af.algorithm.where(x_coords<left_boundary)
  wall_x_right = af.algorithm.where(x_coords>right_boundary)
  wall_y_bot   = af.algorithm.where(y_coords<bottom_boundary)
  wall_y_top   = af.algorithm.where(y_coords>top_boundary)

  n1  = af.Array.elements(wall_x_left)
  n2  = af.Array.elements(wall_x_right)

  x1l = af.randu(n1) 
  x2l = af.randu(n1)
  x3l = af.randu(n1)
  
  x1r = af.randu(n2)
  x2r = af.randu(n2)
  x3r = af.randu(n2)


  sol[wall_x_left]  = 0
  sol[wall_x_right] = 1 - 1e-12
  
  sol[wall_x_left + 2*no_of_particles] = af.arith.sqrt(-2*T_left*af.arith.log(x1l))
  sol[wall_x_right+ 2*no_of_particles] = af.arith.sqrt(-2*T_right*af.arith.log(x1r))*(-1)

  sol[wall_x_left + 3*no_of_particles] = np.sqrt(T_left)*af.arith.sqrt(-2*af.arith.log(x2l))*af.arith.cos(2*np.pi*x3l)
  sol[wall_x_right+ 3*no_of_particles] = np.sqrt(T_right)*af.arith.sqrt(-2*af.arith.log(x2r))*af.arith.cos(2*np.pi*x3r)

  sol[wall_y_bot]   = sol[wall_y_bot] + 1
  sol[wall_y_top]   = sol[wall_y_top] - 1
    
for i in range(no_of_particles):
	j    = sol[(i+1):no_of_particles]        
	k    = sol[(i+1+no_of_particles):2*no_of_particles]
	velx = af.data.constant(af.algorithm.sum(sol[i+2*no_of_particles]),af.Array.elements(j))
	vely = af.data.constant(af.algorithm.sum(sol[i+3*no_of_particles]),af.Array.elements(k))
    
    velx_others = sol[(i+1+2*no_of_particles):3*no_of_particles]
    vely_others = sol[(i+1+3*no_of_particles):4*no_of_particles]
    
    x_particle = af.data.constant(af.algorithm.sum(sol[i]),af.Array.elements(j))
    y_particle = af.data.constant(af.algorithm.sum(sol[i+no_of_particles]),af.Array.elements(k))        
    j          = j-x_particle
    k          = k-y_particle
    dist       = af.arith.sqrt(j**2+k**2)
    j          = j/dist
    k          = k/dist
    
    test_collision = dist<0.01 
    indices        = af.algorithm.where(test_collision)
    
    if(np.sum(test_collision)!=0):
      p           = (velx*j+vely*k-velx_others*j-vely_others*k)*test_collision
      velx        = velx*test_collision-p*j
      velx_others = velx_others+p*j
      vely        = vely*test_collision-p*k
      vely_others = vely_others+p*k
      index       = np.random.randint(0,af.Array.elements(indices))

      sol[i+1+indices[index]+2*no_of_particles] = velx_others[indices[index]]
      sol[i+1+indices[index]+3*no_of_particles] = vely_others[indices[index]]

      sol[i+2*no_of_particles] = velx[indices[index]]
      sol[i+3*no_of_particles] = vely[indices[index]]


  old=sol
  v_x = af.algorithm.sum(sol[2*no_of_particles:3*no_of_particles])
  v_y = af.algorithm.sum(sol[3*no_of_particles:4*no_of_particles])
  En  = af.algorithm.sum(sol[2*no_of_particles:3*no_of_particles]**2 + sol[3*no_of_particles:4*no_of_particles]**2)
  qx  = af.algorithm.sum(sol[2*no_of_particles:3*no_of_particles]*(sol[2*no_of_particles:3*no_of_particles]**2 + sol[3*no_of_particles:4*no_of_particles]**2))

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