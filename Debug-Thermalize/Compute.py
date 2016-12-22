import numpy as np
np.set_printoptions(threshold=np.nan)
import h5py
from SimulationData import *

h5f = h5py.File('Initial.h5', 'r')
initial_conditions = h5f['initial_conditions'][:]
time = h5f['time'][:]
h5f.close()

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

""" Initialization of solution matrix for all particles """

altr=altl=0
pressuredata = np.zeros(time.size,dtype=np.float)
heatfluxdata_x = np.zeros(time.size,dtype=np.float)
heatfluxdata_y = np.zeros(time.size,dtype=np.float)
heatfluxdata_z = np.zeros(time.size,dtype=np.float)
tempassigned=0
counter=0
tempglobal=0

""" Solving """
sol1 = np.zeros(6*no_of_particles,dtype=np.float)
sol2 = np.zeros(6*no_of_particles,dtype=np.float)
sol = np.zeros(6*no_of_particles,dtype=np.float)
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
    sol1=sol   
    
    while(1):       
        for i in range(no_of_particles):
            x_1=np.random.rand(1)
            x_2=np.random.rand(1)
            y_1=const2*np.sqrt(-2*np.log(x_1))*np.cos(2*np.pi*x_2)
            y_2=const2*np.sqrt(-2*np.log(x_1))*np.sin(2*np.pi*x_2)
            x_3=np.random.rand(1)
            x_4=np.random.rand(1)
            y_3=const2*np.sqrt(-2*np.log(x_3))*np.cos(2*np.pi*x_4)
            y_4=const2*np.sqrt(-2*np.log(x_3))*np.sin(2*np.pi*x_4)
            x_5=np.random.rand(1)
            x_6=np.random.rand(1)
            y_5=const2*np.sqrt(-2*np.log(x_5))*np.cos(2*np.pi*x_6)
            y_6=const2*np.sqrt(-2*np.log(x_5))*np.sin(2*np.pi*x_6)
                    
            if(sol[i]>=1):
                sol[i]=1
                
                

                if(altr%2==0):
                    sol[i+3*no_of_particles]=abs(y_1)*(-1)
                    sol[i+4*no_of_particles]=abs(y_2)*np.random.choice([1,-1])
                    sol[i+5*no_of_particles]=abs(y_3)*np.random.choice([1,-1])
                else:        
                    sol[i+3*no_of_particles]=abs(y_4)*(-1)
                    sol[i+4*no_of_particles]=abs(y_5)*np.random.choice([1,-1])
                    sol[i+5*no_of_particles]=abs(y_6)*np.random.choice([1,-1])
                altr=altr+1
                tempassigned=tempassigned+sol[i+3*no_of_particles]**2+sol[i+4*no_of_particles]**2+sol[i+5*no_of_particles]**2
                counter=counter+1
                sol2=sol
                print(sol2-sol1)        
                prompt=input("Hit Enter")

                
            if(sol[i]<=0):
                sol[i]=0
                
                
                if(altl%2==0):
                    sol[i+3*no_of_particles]=abs(y_1)
                    sol[i+4*no_of_particles]=abs(y_2)*np.random.choice([1,-1])
                    sol[i+5*no_of_particles]=abs(y_3)*np.random.choice([1,-1])
                else:        
                    sol[i+3*no_of_particles]=abs(y_4)
                    sol[i+4*no_of_particles]=abs(y_5)*np.random.choice([1,-1])
                    sol[i+5*no_of_particles]=abs(y_6)*np.random.choice([1,-1])
                altl=altl+1
                tempassigned=tempassigned+sol[i+3*no_of_particles]**2+sol[i+4*no_of_particles]**2+sol[i+5*no_of_particles]**2
                counter=counter+1
                sol2=sol
                print(sol2-sol1)        
                prompt=input("Hit Enter")

                
            if(i==(no_of_particles-1)):
                tempassigned=tempassigned/(3*counter)
        sol2=sol
        print(sol2-sol1)        
        print(tempassigned)
        prompt=input("Hit Enter")
        if(abs(tempassigned-T_walls)<0.1):
            break
        altl=alr=tempassigned=counter=0   
    
    sol2=sol

    print(sol2-sol1)

    for i in range(no_of_particles,3*no_of_particles):
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

