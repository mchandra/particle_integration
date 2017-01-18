import numpy as np
import h5py

no_of_particles=100
m=5
left_boundary = 0
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x = length_of_box_x*np.random.rand(no_of_particles)

bottom_boundary = 0
top_boundary = 1
length_of_box_y           = length_of_box_x
initial_conditions_position_y = length_of_box_y*np.random.rand(no_of_particles)

initial_conditions_velocity_x= np.random.rand(no_of_particles)*np.random.choice([-1,1],size=no_of_particles)
initial_conditions_velocity_y= np.random.rand(no_of_particles)*np.random.choice([-1,1],size=no_of_particles)

initial_conditions = np.concatenate([initial_conditions_position_x,initial_conditions_position_y,\
                                    initial_conditions_velocity_x,initial_conditions_velocity_y],axis = 0)

box_crossing_time_scale = length_of_box_x/np.max(initial_conditions_velocity_x)
final_time            = 200
dt   = 0.00005
time = np.arange(0, final_time, dt)
mom_x = np.zeros(time.size)
mom_y = np.zeros(time.size)
Energy = np.zeros(time.size)
KEnergy = np.zeros(time.size)

def potential(a,x):
    potential= 200*(-1*np.tanh(a*x)+1)
    return(potential)

def diffpotential(a,x):
    order1=1e10*(potential(a,x+1e-10)-potential(a,x))
    order2=5e9*(potential(a,x+1e-10)-potential(a,x-1e-10))
    order4=(8*(potential(a,x+1e-10)-potential(a,x-1e-10))+potential(a,x-2e-10)-potential(a,x+2e-10))/(12e-10)
    return(order4)

def calcEnergy(asol):
    TEnergy=0
    for i in range(no_of_particles):
        x=asol[0:no_of_particles]
        y=asol[no_of_particles:2*no_of_particles]
        x=np.delete(x,i)-x[i]
        y=np.delete(y,i)-y[i]
        dist=np.sqrt(x**2+y**2)
        TEnergy+=np.sum(potential(300,dist))
    return(TEnergy)


def Verlet(a,initial_conditions,dt):

    for i in range(no_of_particles):
        x=initial_conditions[0:no_of_particles].copy()
        y=initial_conditions[no_of_particles:2*no_of_particles].copy()
        v_x=initial_conditions[2*no_of_particles:3*no_of_particles].copy()
        v_y=initial_conditions[3*no_of_particles:4*no_of_particles].copy()

        x=np.delete(x,i)-x[i]
        y=np.delete(y,i)-y[i]
        vector = np.array([x,y])
        dist=np.sqrt(x**2+y**2)
        nvector = vector/dist

        F_x = np.sum(diffpotential(a,dist)*nvector[0])
        v_x[i] = v_x[i] + 0.5*(F_x/m)*dt

        F_y = np.sum(diffpotential(a,dist)*nvector[1])
        v_y[i] = v_y[i] + 0.5*(F_y/m)*dt

    x=initial_conditions[0:no_of_particles].copy()
    y=initial_conditions[no_of_particles:2*no_of_particles].copy()

    x_new = x + v_x*dt
    y_new = y + v_y*dt

    for i in range(no_of_particles):
        temp_x=x_new.copy()
        temp_y=y_new.copy()
        x=np.delete(temp_x,i)-temp_x[i]
        y=np.delete(temp_y,i)-temp_y[i]
        vector = np.array([x,y])
        dist=np.sqrt(x**2+y**2)
        nvector = vector/dist

        F_x = np.sum(diffpotential(a,dist)*nvector[0])
        v_x[i] = v_x[i] + 0.5*(F_x/m)*dt

        F_y = np.sum(diffpotential(a,dist)*nvector[1])
        v_y[i] = v_y[i] + 0.5*(F_y/m)*dt

    nextstep=np.concatenate([x_new, y_new, v_x, v_y],axis=0)
    return(nextstep)

sol = np.zeros(4*no_of_particles,dtype=np.float)
old= np.zeros(4*no_of_particles,dtype=np.float)

""" Solving """


for time_index,t0 in enumerate(time):
    t0 = time[time_index]
    print("Computing For Time Index = ",time_index)
    if(time_index==time.size-1):
        break
    t1 = time[time_index+1]
    t = [t0, t1]
    if(time_index==0):
        initial_conditions = initial_conditions
    else:
        initial_conditions = old
    
    sol = Verlet(300,initial_conditions,dt)

    x_coords=sol[0:no_of_particles]
    y_coords=sol[no_of_particles:2*no_of_particles]

    wall_y_bot=np.where(y_coords<bottom_boundary)
    wall_y_top=np.where(y_coords>top_boundary)
    wall_x_left=np.where(x_coords<left_boundary)
    wall_x_right=np.where(x_coords>right_boundary)

    y_coords[wall_y_bot[0]]=y_coords[wall_y_bot[0]]+1
    y_coords[wall_y_top[0]]=y_coords[wall_y_top[0]]-1
    x_coords[wall_x_left[0]]=x_coords[wall_x_left[0]]+1
    x_coords[wall_x_right[0]]=x_coords[wall_x_left[0]]-1

    sol[no_of_particles+wall_y_bot[0]]=y_coords[wall_y_bot[0]]
    sol[no_of_particles+wall_y_top[0]]=y_coords[wall_y_top[0]]
    sol[wall_x_left[0]]=x_coords[wall_x_left[0]]
    sol[wall_x_right[0]]=x_coords[wall_x_right[0]]

    old=sol
    p_x = np.sum(sol[2*no_of_particles:3*no_of_particles])
    p_y = np.sum(sol[3*no_of_particles:4*no_of_particles])
    En = np.sum(sol[2*no_of_particles:3*no_of_particles]**2 + sol[3*no_of_particles:4*no_of_particles]**2)

    mom_x[time_index] = p_x
    mom_y[time_index] = p_y
    KEnergy[time_index] = 0.5*En
    Energy[time_index] = 0.5*En + calcEnergy(sol)
    
    if((time_index%1000)==0):
        h5f = h5py.File('solution_data/solution_'+str(time_index)+'.h5', 'w')
        h5f.create_dataset('sol', data=sol)
        h5f.create_dataset('mom_x', data=mom_x)
        h5f.create_dataset('mom_y', data=mom_y)
        h5f.create_dataset('KEnergy', data=KEnergy)
        h5f.create_dataset('Energy', data=Energy)
        h5f.create_dataset('time', data=time)
        h5f.close()





