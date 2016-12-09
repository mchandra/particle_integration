import numpy as np
import pylab as pl
import h5py

pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
#pl.rcParams['text.usetex']     = True
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'

pl.rcParams['xtick.major.size'] = 8     
pl.rcParams['xtick.minor.size'] = 4     
pl.rcParams['xtick.major.pad']  = 8     
pl.rcParams['xtick.minor.pad']  = 8     
pl.rcParams['xtick.color']      = 'k'     
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'    

pl.rcParams['ytick.major.size'] = 8     
pl.rcParams['ytick.minor.size'] = 4     
pl.rcParams['ytick.major.pad']  = 8     
pl.rcParams['ytick.minor.pad']  = 8     
pl.rcParams['ytick.color']      = 'k'     
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'   


""" Setting number of particles and other parameters"""

no_of_particles = 1000
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
Tc = 50
Th = 50
T = 10
const_right=np.sqrt((k*Tc)/(m))
const_left=np.sqrt((k*Th)/(m))
const = np.sqrt((k*T)/(m))
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
final_time            = 60 * box_crossing_time_scale
dt   = 0.008 * box_crossing_time_scale
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
heatfluxdatax = np.zeros(time.size,dtype=np.float)
heatfluxdatay = np.zeros(time.size,dtype=np.float)
heatfluxdataz = np.zeros(time.size,dtype=np.float)

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
        x_1=np.random.rand(1)
        x_2=np.random.rand(1)
        x_3=np.random.rand(1)
        x_4=np.random.rand(1)
        x_5=np.random.rand(1)
        x_6=np.random.rand(1)
        if(sol[i]>=right_boundary):        #Cold Resevoir right
            if(alternator%2==0):
                sol[3*no_of_particles+i] = abs(np.sqrt(const_right)*np.sqrt(-2*np.log(x_1))*np.cos(2*np.pi*x_2)) * (-1)
                sol[4*no_of_particles+i] = abs(np.sqrt(const_right)*np.sqrt(-2*np.log(x_3))*np.cos(2*np.pi*x_4))
                sol[5*no_of_particles+i] = abs(np.sqrt(const_right)*np.sqrt(-2*np.log(x_5))*np.cos(2*np.pi*x_6))
            else:
                sol[3*no_of_particles+i] = abs(np.sqrt(const_right)*np.sqrt(-2*np.log(x_1))*np.sin(2*np.pi*x_2)) * (-1)
                sol[4*no_of_particles+i] = abs(np.sqrt(const_right)*np.sqrt(-2*np.log(x_3))*np.sin(2*np.pi*x_4))
                sol[5*no_of_particles+i] = abs(np.sqrt(const_right)*np.sqrt(-2*np.log(x_5))*np.sin(2*np.pi*x_6))
            alternator=alternator+1
        if(sol[i]<=left_boundary):       #Hot Resevoir left side
            if(alternatol%2==0):
                sol[3*no_of_particles+i] = abs(np.sqrt(const_left)*np.sqrt(-2*np.log(x_1))*np.cos(2*np.pi*x_2)) * (+1)
                sol[4*no_of_particles+i] = abs(np.sqrt(const_left)*np.sqrt(-2*np.log(x_3))*np.cos(2*np.pi*x_4))
                sol[5*no_of_particles+i] = abs(np.sqrt(const_left)*np.sqrt(-2*np.log(x_5))*np.cos(2*np.pi*x_6))
            else:
                sol[3*no_of_particles+i] = abs(np.sqrt(const_left)*np.sqrt(-2*np.log(x_1))*np.sin(2*np.pi*x_2)) * (+1)
                sol[4*no_of_particles+i] = abs(np.sqrt(const_left)*np.sqrt(-2*np.log(x_3))*np.sin(2*np.pi*x_4))
                sol[5*no_of_particles+i] = abs(np.sqrt(const_left)*np.sqrt(-2*np.log(x_5))*np.sin(2*np.pi*x_6))
            alternatol=alternatol+1
   
    for i in range(2*no_of_particles,3*no_of_particles):
        if(sol[i]>=right_boundary):
            sol[i] = sol[i] - length_of_box_x
        if(sol[i]<=left_boundary):
            sol[i] = sol[i] + length_of_box_x
    for i in range(no_of_particles):
        pressuredata[time_index] +=  sol[3*no_of_particles+i]**2+sol[4*no_of_particles+i]**2 +sol[5*no_of_particles+i]**2
        heatfluxdatax[time_index] += (sol[3*no_of_particles+i]**2+sol[4*no_of_particles+i]**2 +sol[5*no_of_particles+i]**2)*sol[3*no_of_particles+i]
        heatfluxdatay[time_index] += (sol[3*no_of_particles+i]**2+sol[4*no_of_particles+i]**2 +sol[5*no_of_particles+i]**2)*sol[4*no_of_particles+i]
        heatfluxdataz[time_index] += (sol[3*no_of_particles+i]**2+sol[4*no_of_particles+i]**2 +sol[5*no_of_particles+i]**2)*sol[5*no_of_particles+i]
    old=sol
    #h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'w')
    #h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=sol)
    #h5f.close()
    #print(sol[0])

for i in range(time.size):
    pressuredata[i] = pressuredata[i]/(3*no_of_particles)
    heatfluxdatax[i] = heatfluxdatax[i]/(3*no_of_particles)
    heatfluxdatay[i] = heatfluxdatay[i]/(3*no_of_particles)
    heatfluxdataz[i] = heatfluxdataz[i]/(3*no_of_particles)
pl.figure()
pl.plot(time,heatfluxdatax,'r--',label='heat flux in x')
#pl.legend()
pl.plot(time,heatfluxdatay,'g--',label='heat flux in y')
#pl.legend()
pl.plot(time,heatfluxdataz,'b--',label='heat flux in z')
#pl.legend()
pl.plot(time,pressuredata,'y',label='sigma v^2')
#pl.legend()
pl.savefig('time_plot.png')
pl.xlabel('time')
pl.ylabel('pressure and heat flux')
pl.legend()
pl.clf()
