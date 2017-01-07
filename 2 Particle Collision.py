import numpy as np
import pylab as pl
import h5py

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
y_1=np.linspace(50,55,101)
angle1=np.zeros(101)
angle2=np.zeros(101)
impactf=np.linspace(0,1,101)

x_1=43
y_1=50
m_1=5
r_1=5
v_x_1=5
v_y_1=0.01

x_2=50
y_2=50
m_2=300
r_2=300
v_x_2=0
v_y_2=0

m=[m_1,-1*m_2] #To make the forces on the particles, act in opposite directions
initial_conditions=np.array([x_1,x_2,y_1,y_2,v_x_1,v_x_2,v_y_1,v_y_2])

def V(x):
	V=200*(-1*np.tanh(5*(x-5.2))+1)
	return(V)

def Verlet(initial_conditions,t):

	x=initial_conditions[0:no_of_particles]
	y=initial_conditions[no_of_particles:2*no_of_particles]
	v_x=initial_conditions[2*no_of_particles:3*no_of_particles]
	v_y=initial_conditions[3*no_of_particles:4*no_of_particles]
	x_new = x + v_x*(t[1]-t[0])
	y_new = y + v_y*(t[1]-t[0])

	vector = np.array([x_new[1]-x_new[0],y_new[1]-y_new[0]])
	nvector = vector/np.linalg.norm(vector)

	F_x = 500000*nvector[0]*(V(np.linalg.norm(vector+0.000001)) - V(np.linalg.norm(vector-0.000001)))
	v_x_new = v_x + (F_x/m)*(t[1]-t[0])

	F_y = 500000*nvector[1]*(V(np.linalg.norm(vector+0.000001)) - V(np.linalg.norm(vector-0.000001)))
	v_y_new = v_y + (F_y/m)*(t[1]-t[0])
	nextstep=np.concatenate([x_new, y_new, v_x_new, v_y_new],axis=0)
	return(nextstep)

count=0



""" Discretizing time and making sure scaling is done right """

final_time            = 1
dt   = 0.000005
time = np.arange(0, final_time, dt)
energydata = np.zeros(time.size,dtype=np.float) 
kedata_1 = np.zeros(time.size,dtype=np.float) 
kedata_2 = np.zeros(time.size,dtype=np.float) 
pxdata = np.zeros(time.size,dtype=np.float) 
pydata = np.zeros(time.size,dtype=np.float)

""" Initialization of solution matrix for all particles """

sol = np.zeros(4*no_of_particles,dtype=np.float)
fsol = np.zeros(4*no_of_particles,dtype=np.float)
""" Solving """

old= np.zeros(4*no_of_particles,dtype=np.float)
for time_index,t0 in enumerate(time):
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

	
	old=sol
	momentum_x=m_1*sol[4]+m_2*sol[5]	
	momentum_y=m_1*sol[6]+m_2*sol[7]
	energy=0.5*m_1*(sol[4]**2+sol[6]**2)+0.5*m_2*(sol[5]**2+sol[7]**2)+V(np.sqrt((sol[3]-sol[2])**2+(sol[1]-sol[0])**2))
	ke_1=0.5*m_1*(sol[4]**2+sol[6]**2)
	ke_2=0.5*m_2*(sol[5]**2+sol[7]**2)
	energydata[time_index]=energy
	kedata_1[time_index]=ke_1
	kedata_2[time_index]=ke_2
	pxdata[time_index]=momentum_x
	pydata[time_index]=momentum_y


pl.plot(time,energydata,label='$\mathrm{Total}$ $\mathrm{Energy}$')
#pl.plot(time,kedata_1,label='$\mathrm{Kinetic}$ $\mathrm{Energy}$ $\mathrm{of}$ $\mathrm{Particle1}$')
#pl.plot(time,kedata_2,label='$\mathrm{Kinetic}$ $\mathrm{Energy}$ $\mathrm{of}$ $\mathrm{Particle2}$')
pl.plot(time,pxdata,label='$\mathrm{Momentum}$ $\mathrm{in}$ $x$')
pl.plot(time,pydata,label='$\mathrm{Momentum}$ $\mathrm{in}$ $y$')


pl.title('$\mathrm{Energy/Momentum}$ $\mathrm{Variation}$')
pl.xlabel('$\mathrm{Time}$')
pl.ylabel('$\mathrm{Energy/Momentum}$')
pl.legend(loc='center right')
pl.savefig('Variation2Particle.png') 
