import numpy as np
import h5py
from scipy.special import erfinv
import numexpr as ne

""" Setting number of particles and other parameters"""
no_of_particles = 1000000
x_div = 100
y_div = 100


""" Initial conditions """

left_boundary                 = 0
right_boundary                = 1.0
length_of_box_x               = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x*np.random.rand(no_of_particles)

bottom_boundary               = 0
top_boundary                  = 1.0
length_of_box_y               = top_boundary - bottom_boundary
initial_conditions_position_y = left_boundary + length_of_box_x*np.random.rand(no_of_particles)

dx = (length_of_box_x/x_div)
x  = np.arange(left_boundary, right_boundary,dx)
x  = np.concatenate((x,[right_boundary]),axis = 0)

dy = (length_of_box_y/y_div)
y  = np.arange(left_boundary, right_boundary,dy)
y  = np.concatenate((y,[top_boundary]),axis = 0)

""" Setting velocities according to maxwellian distribution """

k=1.0
m=1.0
T=1.5
T_walls=2.0

initial_conditions_velocity_x=np.zeros(no_of_particles,dtype=np.float)
initial_conditions_velocity_y=np.zeros(no_of_particles,dtype=np.float)

x1 = np.random.rand(no_of_particles)
x2 = np.random.rand(no_of_particles)


initial_conditions_velocity_x=np.sqrt(2*T)*erfinv(2*x1-1)
initial_conditions_velocity_y=np.sqrt(2*T)*erfinv(2*x2-1)

""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                  initial_conditions_velocity_x,initial_conditions_velocity_y],axis = 0)

""" Discretizing time and making sure scaling is done right """

box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)
final_time              = 5.0
dt   = 0.05 * box_crossing_time_scale
time = np.arange(0, final_time, dt)
print(time.size)
temp_global = np.zeros(100)


"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,t):

  x=initial_conditions[0:no_of_particles]
  y=initial_conditions[no_of_particles:2*no_of_particles]
  v_x=initial_conditions[2*no_of_particles:3*no_of_particles]
  v_y=initial_conditions[3*no_of_particles:4*no_of_particles]
  dt=t[1]-t[0]
  x_new = x + v_x*dt
  v_x_new = v_x
  y_new = y + v_y*dt
  v_y_new = v_y
  nextstep=np.concatenate([x_new, y_new, v_x_new, v_y_new],axis=0)
  return(nextstep)

  return(sol) 
""" Initialization of solution matrix for all particles """

sol = np.zeros(4*no_of_particles,dtype=np.float)

pressuredata = np.zeros(time.size,dtype=np.float)
heatfluxdata_x = np.zeros(time.size,dtype=np.float)
heatfluxdata_y = np.zeros(time.size,dtype=np.float)

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
  v_x=sol[2*no_of_particles:3*no_of_particles]
  v_y=sol[3*no_of_particles:4*no_of_particles]

  """Implementing the Thermal B.C's"""

  collided_right=np.where(x_coords>1)
  collided_left=np.where(x_coords<0)

  x1 = np.random.rand(collided_right[0].size)
  x2 = np.random.rand(collided_right[0].size)

  x_coords[collided_right[0]] = 1 - 1e-10
  v_x[collided_right[0]] = np.sqrt(-2*T_walls*np.log(x1))*(-1)    
  v_y[collided_right[0]] = np.sqrt(2*T_walls)*erfinv(2*x2-1)
  
  sol[collided_right[0]+2*no_of_particles] = v_x[collided_right[0]]
  sol[collided_right[0]+3*no_of_particles] = v_y[collided_right[0]]
  sol[collided_right[0]] = x_coords[collided_right[0]]
  
  x1 = np.random.rand(collided_left[0].size)
  x2 = np.random.rand(collided_left[0].size)
   
  x_coords[collided_left[0]] = 0
  v_x[collided_left[0]] = np.sqrt(-2*T_walls*np.log(x1))  
  v_y[collided_left[0]] = np.sqrt(2*T_walls)*erfinv(2*x2-1)   

  sol[2*no_of_particles+collided_left[0]] = v_x[collided_left[0]]
  sol[3*no_of_particles+collided_left[0]] = v_y[collided_left[0]]
  sol[collided_left[0]] = x_coords[collided_left[0]]

  """Implementation of Periodic B.C's"""

  wall_y_bot=np.where(y_coords<bottom_boundary)
  wall_y_top=np.where(y_coords>top_boundary)
  
  y_coords[wall_y_bot[0]]=y_coords[wall_y_bot[0]]+1
  y_coords[wall_y_top[0]]=y_coords[wall_y_top[0]]-1
  
  sol[no_of_particles+wall_y_bot[0]]=y_coords[wall_y_bot[0]]
  sol[no_of_particles+wall_y_top[0]]=y_coords[wall_y_top[0]]

  x_zones   = 100 * sol[0:no_of_particles]
  y_zones   = 100 * sol[no_of_particles:2*no_of_particles]
  x_zones   = x_zones.astype(int)
  y_zones   = y_zones.astype(int)
  zone      = 100*y_zones + x_zones
  zonecount = np.bincount(zone)
  
  temp = np.zeros(zonecount.size)

  for i in range(10000):
    indices = np.where(zone == i)[0]
    temp[i] = 0.5*np.sum(sol[indices+2*no_of_particles]**2 + sol[indices+3*no_of_particles]**2)
  
  temp=temp/zonecount
    
  for i in range(10000):
    indices = np.where(zone == i)[0]
    x1 = np.random.rand(zonecount[i])
    x2 = np.random.rand(zonecount[i])
    sol[indices+2*no_of_particles] = np.sqrt(2*temp[i])*erfinv(2*x1-1) 
    sol[indices+3*no_of_particles] = np.sqrt(2*temp[i])*erfinv(2*x2-1)

  temp_global = np.sum(temp, axis=0)/100
  old=sol

  pressure=0
  heatflux_x=0
  heatflux_y=0

  pressure=pressure+ne.evaluate("sum(v_x**2+v_y**2)")
  heatflux_x=heatflux_x+ne.evaluate("sum(v_x*(v_x**2+v_y**2))")
  heatflux_y=heatflux_y+ne.evaluate("sum(v_y*(v_x**2+v_y**2))")

  heatflux_x=heatflux_x/no_of_particles
  heatflux_y=heatflux_y/no_of_particles
  pressure=pressure/no_of_particles

  pressuredata[time_index]=pressure
  heatfluxdata_x[time_index]=heatflux_x
  heatfluxdata_y[time_index]=heatflux_y


  filename = 'Files/Collisional/post' + '%04d'%time_index + '.h5'
  h5f = h5py.File(filename, 'w')
  h5f.create_dataset('time', data=time)
  h5f.create_dataset('temp', data=temp_global)
  h5f.create_dataset('sol', data=sol)
  h5f.create_dataset('heatflux_x', data=heatfluxdata_x)
  h5f.create_dataset('heatflux_y', data=heatfluxdata_y)
  h5f.create_dataset('pressure', data=pressuredata)
  h5f.close()