import numpy as np
import h5py
import pylab as pl
import sympy
pl.axes().set_aspect('equal', 'datalim')

pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 15  
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

m=[5 ,-300]

def potential(a,x):
	potential=200*(-1*np.tanh(a*(x-5))+1)
	return(potential)

def diffpotential(a,x):
	order1=1e10*(potential(a,x+1e-10)-potential(a,x))
	order2=5e9*(potential(a,x+1e-10)-potential(a,x-1e-10))
	order4=(8*(potential(a,x+1e-10)-potential(a,x-1e-10))+potential(a,x-2e-10)-potential(a,x+2e-10))/(12e-10)
	return(order4)


def Verlet(a,initial_conditions,dt):

	x=initial_conditions[0:no_of_particles]
	y=initial_conditions[no_of_particles:2*no_of_particles]
	v_x=initial_conditions[2*no_of_particles:3*no_of_particles]
	v_y=initial_conditions[3*no_of_particles:4*no_of_particles]
	vector = np.array([x[1]-x[0],y[1]-y[0]])
	nvector = vector/np.linalg.norm(vector)

	F_x = diffpotential(a,np.linalg.norm(vector))*nvector[0]
	v_x = v_x + 0.5*(F_x/m)*dt

	F_y = diffpotential(a,np.linalg.norm(vector))*nvector[1]
	v_y = v_y + 0.5*(F_y/m)*dt

	x_new = x + v_x*dt
	y_new = y + v_y*dt

	vector = np.array([x_new[1]-x_new[0],y_new[1]-y_new[0]])
	nvector = vector/np.linalg.norm(vector)

	F_x = diffpotential(a,np.linalg.norm(vector))*nvector[0]
	v_x = v_x + 0.5*(F_x/m)*dt

	F_y = diffpotential(a,np.linalg.norm(vector))*nvector[1]
	v_y = v_y + 0.5*(F_y/m)*dt

	nextstep=np.concatenate([x_new, y_new, v_x, v_y],axis=0)
	return(nextstep)

no_of_particles=2
impactf=np.linspace(50,55,101)
angle1 = np.zeros(101,dtype=np.float)
count = 0

x_1=44.5
y_1=54
m_1=5
r_1=5
v_x_1=1
v_y_1=0

x_2=50
y_2=50
m_2=300
r_2=300
v_x_2=0
v_y_2=0

initial_conditions=np.array([x_1,x_2,y_1,y_2,v_x_1,v_x_2,v_y_1,v_y_2])

""" Discretizing time and making sure scaling is done right """

final_time            = 6
time=np.zeros(1)

x_stored=x_1
y_stored=y_1
vx_stored=v_x_1
vy_stored=v_y_1
error=0

""" Initialization of solution matrix for all particles """

sol = np.zeros(4*no_of_particles,dtype=np.float)

""" Solving """
time_index=0
old= np.zeros(4*no_of_particles,dtype=np.float)
while(time[-1]<=6):
	t0 = time[time_index]
	print("Computing for Time Index = ",time_index," Time = ",t0)
	testdt=0.01
	if(time_index==0):
		initial_conditions = initial_conditions
	else:
		initial_conditions = old

	while(1):
		sol1 = Verlet(300,initial_conditions,testdt)
		testdt=testdt/2
		sol2 = Verlet(300,initial_conditions,testdt)
		errorstep= abs(np.max(sol1-sol2))
		if(errorstep<1e-5):
			error=np.append(error,[errorstep])
			break
	old=sol1
	t1 = t0 + testdt
	time = np.append(time,[t1])

	x_stored = np.append(x_stored,[sol[0]])
	y_stored = np.append(y_stored,[sol[2]])
	vx_stored = np.append(vx_stored,[sol[4]])
	vy_stored = np.append(vy_stored,[sol[6]])
	time_index = time_index + 1

if(sol[4]<0):
	angle1[count] = np.pi - np.arctan(sol[6]/abs(sol[4]))
else:
	angle1[count] = np.arctan(sol[6]/sol[4])
#count=count+1

h5f = h5py.File('potential2.h5', 'w')
h5f.create_dataset('x', data=x_stored)
h5f.create_dataset('y', data=y_stored)
h5f.create_dataset('vx', data=vx_stored)
h5f.create_dataset('vy', data=vy_stored)
h5f.create_dataset('time', data=time)
h5f.create_dataset('error', data=error)
h5f.close()

h5f = h5py.File('v1-2.h5', 'w')
h5f.create_dataset('angle1', data=angle1)
h5f.close()