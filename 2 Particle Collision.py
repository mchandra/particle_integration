import numpy as np
import pylab as pl

pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
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

no_of_particles=2

x_1=20
y_1=50
m_1=5
r_1=5
v_x_1=200
v_y_1=0

x_2=50
y_2=50
m_2=300
r_2=300
v_x_2=0
v_y_2=0

initial_conditions=np.array([x_1,x_2,y_1,y_2,v_x_1,v_x_2,v_y_1,v_y_2])

""" Discretizing time and making sure scaling is done right """

final_time            = 0.25
dt   = 0.0005
time = np.arange(0, final_time, dt)

"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,t):
    x=initial_conditions[0:no_of_particles]
    y=initial_conditions[no_of_particles:2*no_of_particles]
    v_x=initial_conditions[2*no_of_particles:3*no_of_particles]
    v_y=initial_conditions[3*no_of_particles:4*no_of_particles]
    x_new = x + v_x*(t[1]-t[0])
    v_x_new = v_x 
    y_new = y + v_y*(t[1]-t[0])
    v_y_new = v_y 
    nextstep=np.concatenate([x_new, y_new, v_x_new, v_y_new],axis=0)
    return(nextstep)

""" Initialization of solution matrix for all particles """

sol = np.zeros(4*no_of_particles,dtype=np.float)

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
        
    for i in range(0,2*no_of_particles):
        if(sol[i]>=99 or sol[i]<=1):
            sol[i+2*no_of_particles] = (-1)*sol[i+2*no_of_particles] 

    """Accounting of collisions between the two particles"""
    old=sol
    dist=np.sqrt((sol[1]-sol[0])**2+(sol[3]-sol[2])**2)
    print(dist)

    if (dist<5.0):
        normal=np.array([(sol[1]-sol[0]),(sol[3]-sol[2])])
        normal=normal/dist
        p=2*(sol[4]*normal[0]+sol[6]*normal[1]-sol[5]*normal[0]-sol[7]*normal[1])/(m_1+m_2)
        sol[4]=sol[4]-p*m_2*normal[0]
        sol[5]=sol[5]+p*m_1*normal[0]
        sol[6]=sol[6]-p*m_2*normal[1]
        sol[7]=sol[7]+p*m_1*normal[1]
        print(sol[4],sol[5],sol[6],sol[7])

    pl.plot(sol[0],sol[2], 'o',color='blue', markersize=5, alpha = 1)
    pl.plot(sol[1],sol[3], 'o',color='blue', markersize=60, alpha = 0.5)
    pl.xlim(0, 100)
    pl.ylim(0, 100)
    pl.title('$\mathrm{Time}$ = ' + str(time[time_index]) )
    pl.xlabel('$x$')
    pl.ylabel('$y$')
    pl.axes().set_aspect('equal', 'datalim')
    pl.savefig('images/' + '%04d'%time_index + '.png')
    pl.clf() 
