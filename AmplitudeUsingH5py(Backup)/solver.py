from initialization import *


""" Initialization of solution matrix for all particles """


solution_all_current = np.zeros((4*no_of_particles),dtype=np.float)

#%%


""" Solving """
BC='2'


for time_index,t0 in enumerate(time):
    
    t0 = time[time_index] 
    print("Computing for TimeIndex = ",time_index)
    if(time_index==time.size-1):
        break
    t1 = time[time_index+1]
    t = [t0, t1] 
    if(time_index==0):
        h5f = h5py.File('initial_conditions.h5', 'r')
        initial_conditions = h5f['initial_conditions_dataset'][:]
        h5f.close()
       
    else:
        initial_conditions = (solution_all_current[:])
    sol = Verlet(initial_conditions,t)
    solution_all_current[:] = sol
    if(BC=='1'):
    
        for i in range(2*no_of_particles):
            if(solution_all_current[i]>right_boundary):
                solution_all_current[2*no_of_particles+i] = solution_all_current[2*no_of_particles+i] * (-1)
            if(solution_all_current[i]<left_boundary):
                solution_all_current[2*no_of_particles+i] = solution_all_current[2*no_of_particles+i] * (-1)
                
    if(BC=='2'):
        for i in range(2*no_of_particles):
            if(solution_all_current[i]>right_boundary):
                solution_all_current[i] = solution_all_current[i] - length_of_box_x
            if(solution_all_current[i]<left_boundary):
                solution_all_current[i] = solution_all_current[i] + length_of_box_x
                




    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=solution_all_current)
    h5f.close()

