import numpy as np
import pylab as pl


#%%
""" Set plot parameters to make beautiful plots """
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
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

"""Definition of Verlet Integrator"""

def Verlet(initial_conditions,t):
	x=initial_conditions[0]
	v_x=initial_conditions[1]
	x_new = x + v_x*(t[1]-t[0])
	v_x_new = v_x + -1*x_new*(t[1]-t[0]) 
	return([x_new,v_x_new])
	
#%%


""" Initial conditions """

left_boundary = -1
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x*np.random.rand(1)


""" Setting the  velocity distribution"""

initial_conditions_velocity_x =  np.random.uniform(2,5)

""" Combining the initial conditions into one vector"""

initial_conditions = [initial_conditions_position_x,initial_conditions_velocity_x]



""" Discretizing time and making sure scaling is done right """
box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)
final_time            = 10 * box_crossing_time_scale
dt   = 0.001 * box_crossing_time_scale
time = np.arange(0, final_time, dt)

""" Initialization of solution matrix for all particles """
solution_all = np.zeros((2,time.size),dtype=np.float)

""" Solving """
BC=input("1 for Hard,2 for Periodic and any other number for no BC's:")

for time_index,t0 in enumerate(time):
    
    t0 = time[time_index] 
    
    if(time_index==time.size-1):
        break
    t1 = time[time_index+1]
    t = [t0, t1] 
    if(time_index==0):
        initial_conditions = initial_conditions
       
    else:
        initial_conditions = (solution_all[:,time_index-1])
    sol = Verlet(initial_conditions,t)
    solution_all[:,time_index] = sol
    if(BC=='1'):
    
            if(solution_all[0,time_index]>right_boundary):
                solution_all[1,time_index] = solution_all[1,time_index] * (-1)
            if(solution_all[0,time_index]<left_boundary):
                solution_all[1,time_index] = solution_all[1,time_index] * (-1)
                
    if(BC=='2'):
       
            if(solution_all[0,time_index]>right_boundary):
                solution_all[0,time_index] = solution_all[0,time_index] - length_of_box_x
            if(solution_all[0,time_index]<left_boundary):
                solution_all[0,time_index] = solution_all[0,time_index] + length_of_box_x

    else:
        continue


pl.plot(time,solution_all[0,:])
pl.title('Trajectory of Particle')
pl.xlabel('$\mathrm{Time}$')
pl.ylabel('$x$')
pl.show()





