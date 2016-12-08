# Heat flux and temperature variation with x 

from joblib import Parallel, delayed
import multiprocessing
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


""" some initial terms """

no_of_particles = 20000
x_divisions=32
y_divisions=1
z_divisions=1
left_boundary = 0
right_boundary = 1
length_of_box_x = right_boundary - left_boundary


N=5000

num_cores = multiprocessing.cpu_count()
inputs = range(N)   



""" Initializing """

dx = (length_of_box_x/x_divisions)
x = np.arange(left_boundary,right_boundary,dx)
x= np.concatenate((x,[right_boundary]),axis = 0)



def func(time_index):
    print("time index ", str(time_index))
    v_x = np.zeros(no_of_particles)
    v_y = np.zeros(no_of_particles)
    v_z = np.zeros(no_of_particles)
    particle_x_zone = np.zeros((no_of_particles),dtype=np.float)
    v = np.zeros(no_of_particles)
    sigma_v2 = np.zeros(x.size-1)
    sigma_v2_vx = np.zeros(x.size-1)
    sigma_v2_vy = np.zeros(x.size-1)
    sigma_v2_vz = np.zeros(x.size-1)
    count = np.zeros(x.size-1)
    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'r')
    solution_all_current = h5f['solution_all/solution_dataset_'+str(time_index)][:]
    h5f.close()
    v_x = solution_all_current[3*no_of_particles:4*no_of_particles]
    v_y = solution_all_current[4*no_of_particles:5*no_of_particles]
    v_z = solution_all_current[5*no_of_particles:6*no_of_particles]
    for i in range(no_of_particles):
        v[i] = (v_x[i]*v_x[i])+(v_y[i]*v_y[i])+(v_z[i]*v_z[i])
		
    for p in range(0,no_of_particles):
        
        for i in range(x.size-1):
		
            if((solution_all_current[p]>x[i])and(solution_all_current[p]<x[i+1])):
                count[i] += 1
                sigma_v2[i] += v[p]**2
                sigma_v2_vx[i] += v[p]**2*v_x[p]    
                sigma_v2_vy[i] += v[p]**2*v_y[p] 
                sigma_v2_vz[i] += v[p]**2*v_z[p]   
                break
    for i in range(x.size-1):
	    sigma_v2[i] = sigma_v2[i]/count[i]
	    sigma_v2_vx[i] = sigma_v2_vx[i]/count[i]	
	    sigma_v2_vy[i] = sigma_v2_vy[i]/count[i]
	    sigma_v2_vz[i] = sigma_v2_vz[i]/count[i]
    #print(sigma_v2)
    #print(sigma_v2_vx)
    ##print(sigma_v2_vy)
    #print(sigma_v2_vz)
    pl.figure()
    pl.plot(sigma_v2,'r',label='Temperature')
    pl.plot(sigma_v2_vx,'b--',label='heat flux in x')
    pl.plot(sigma_v2_vy,'g--',label='heat flux in y')
    pl.plot(sigma_v2_vz,'y--',label='heat flux in z')
    pl.ylim(-150,200)
    pl.savefig('v2/%04d'%time_index + '.png')
    pl.legend()
    pl.clf()
	



				



if __name__ == '__main__':
    results = Parallel(n_jobs=2)(delayed(func)(i) for i in inputs)




