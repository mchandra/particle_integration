from params import *
from scipy.integrate import odeint


#Equations are:
# dBz/dt = - dEy/dx
# dEy/dt = - dBz/dx


""" Setting number of particles and other parameters"""

no_of_particles = 1

Nx = 300
Ny = 300

dx = Lx/Nx
dy = Ly/Ny

""" Setting the grids """

x_center = np.linspace(-ghost_cells*dx, Lx + ghost_cells*dx, Nx + 1 + 2 * ghost_cells, endpoint=True)
y_center = np.linspace(-ghost_cells*dy, Ly + ghost_cells*dy, Ny + 1 + 2 * ghost_cells, endpoint=True)

""" Setting the offset spatial grids """


x_right = np.linspace(-ghost_cells * dx / 2, Lx + (2 * ghost_cells + 1) * dx / 2, Nx + 1 + 2 * ghost_cells,\
                        endpoint=True\
                     )


y_top = np.linspace(-ghost_cells * dy / 2, Ly + (2 * ghost_cells + 1) * dy / 2, Ny + 1 + 2 * ghost_cells,\
                      endpoint=True\
                   )


""" Initial conditions """

left_boundary = 0
right_boundary = Lx
length_of_box_x = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x * np.random.rand(no_of_particles)
#initial_conditions_position_x = np.ones((no_of_particles), dtype = np.float)*Lx/2

bottom_boundary = 0
top_boundary = Ly
length_of_box_y = length_of_box_x
#initial_conditions_position_y = bottom_boundary + length_of_box_y * np.random.rand(no_of_particles)
initial_conditions_position_y = np.ones((no_of_particles), dtype = np.float)*Ly/2

back_boundary = 0
front_boundary = Ly
length_of_box_z = length_of_box_x
initial_conditions_position_z = back_boundary + length_of_box_z * np.random.rand(no_of_particles)

""" Setting velocities according to maxwellian distribution """

initial_conditions_velocity_x = np.zeros(no_of_particles, dtype=np.float)
initial_conditions_velocity_y = np.zeros(no_of_particles, dtype=np.float)
initial_conditions_velocity_z = np.zeros(no_of_particles, dtype=np.float)

initial_conditions_velocity_y[:] = 2
""" Combining the initial conditions into one vector"""

initial_conditions = np.concatenate([initial_conditions_position_x, initial_conditions_position_y,\
                                     initial_conditions_position_z, initial_conditions_velocity_x,\
                                     initial_conditions_velocity_y, initial_conditions_velocity_z], axis=0)


""" Electric and Magnetic field """

Ez = np.zeros((len(x_center), len(y_center)), dtype=np.float)
Bx = np.zeros((len(x_center), len(y_top)   ), dtype=np.float)
By = np.zeros((len(x_right), len(y_center) ), dtype=np.float)

Bz = np.zeros((len(x_center), len(y_center)), dtype=np.float)
Ex = np.zeros((len(x_center), len(y_top)   ), dtype=np.float)
Ey = np.zeros((len(x_right), len(y_center) ), dtype=np.float)
#print('Bz size is ', Bz.shape)

# mass and charge of particle:

me = 1
charge = 1


""" Writing the spatial grids as a two dimension matrix for vectorization purposes """

X_center_physical, Y_center_physical = np.meshgrid( x_center[ghost_cells:-ghost_cells],\
                                                    y_center[ghost_cells:-ghost_cells]\
                                                  )

""" Writing the offset spatial grids and indices as a two dimension matrix for vectorization purposes """

X_right_physical, Y_top_physical = np.meshgrid( x_right[ghost_cells:-ghost_cells],\
                                                y_top[ghost_cells:-ghost_cells]\
                                              )

I, J = np.meshgrid( range(ghost_cells, len(x_center) - ghost_cells),\
                    range(ghost_cells, len(y_center) - ghost_cells)\
                  )
# Initializing the non relevant fields:

Ey[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = np.sin(2*np.pi*(-X_right_physical))
Bz[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = np.sin(2*np.pi*(-X_right_physical))

#Bz[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = 20

""" Discretizing time and making sure scaling is done right """

# box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)

final_time = 5
dt = np.float(dx / (2 * c))
time = np.arange(0, final_time, dt)

""" Boris method modified Verlet Integrator """


def mag_Verlet(initial_conditions, t, F_interpolated):
    x = initial_conditions[0:no_of_particles]
    y = initial_conditions[no_of_particles:2 * no_of_particles]
    z = initial_conditions[2 * no_of_particles:3 * no_of_particles]
    v_x = initial_conditions[3 * no_of_particles:4 * no_of_particles]
    v_y = initial_conditions[4 * no_of_particles:5 * no_of_particles]
    v_z = initial_conditions[5 * no_of_particles:6 * no_of_particles]
    
    Ex = F_interpolated[0,0:no_of_particles]
    Ey = F_interpolated[1,0:no_of_particles]
    Ez = F_interpolated[2,0:no_of_particles]
    Bx = F_interpolated[3,0:no_of_particles]
    By = F_interpolated[4,0:no_of_particles]
    Bz = F_interpolated[5,0:no_of_particles]

    print('Expected magnetic field = ', np.sin(2*np.pi*(t[0] - x[0])))
    print('Numerical magnetic field = ', Bz )



    x_new = x + v_x * (t[1] - t[0])
    y_new = y + v_y * (t[1] - t[0])
    z_new = z + v_z * (t[1] - t[0])

    v_x_minus = v_x + (charge * Ex * (t[1] - t[0])) / (2 * me)
    v_y_minus = v_y + (charge * Ey * (t[1] - t[0])) / (2 * me)
    v_z_minus = v_z + (charge * Ez * (t[1] - t[0])) / (2 * me)

    t_magx = (charge * Bx * (t[1] - t[0])) / (2 * me)
    t_magy = (charge * By * (t[1] - t[0])) / (2 * me)
    t_magz = (charge * Bz * (t[1] - t[0])) / (2 * me)

    vminus_cross_t_x = (v_y_minus * t_magz) - (v_z_minus * t_magy)
    vminus_cross_t_y = -((v_x_minus * t_magz) - (v_z_minus * t_magx))
    vminus_cross_t_z = (v_x_minus * t_magy) - (v_y_minus * t_magx)

    v_dashx = v_x_minus + vminus_cross_t_x
    v_dashy = v_y_minus + vminus_cross_t_y
    v_dashz = v_z_minus + vminus_cross_t_z

    t_mag = np.sqrt(t_magx ** 2 + t_magy ** 2 + t_magz ** 2)

    s_x = (2 * t_magx) / (1 + abs(t_mag ** 2))
    s_y = (2 * t_magy) / (1 + abs(t_mag ** 2))
    s_z = (2 * t_magz) / (1 + abs(t_mag ** 2))

    v_x_plus = v_x_minus + ((v_dashy * s_z) - (v_dashz * s_y))
    v_y_plus = v_y_minus + -((v_dashx * s_z) - (v_dashz * s_x))
    v_z_plus = v_z_minus + ((v_dashx * s_y) - (v_dashy * s_x))

    v_x_new = v_x_plus + (charge * Ex * (t[1] - t[0])) / (2 * me)
    v_y_new = v_y_plus + (charge * Ey * (t[1] - t[0])) / (2 * me)
    v_z_new = v_z_plus + (charge * Ez * (t[1] - t[0])) / (2 * me)

    nextstep = np.concatenate([x_new, y_new, z_new, v_x_new, v_y_new, v_z_new], axis=0)
    return (nextstep)


#mag_Verlet = np.vectorize(mag_Verlet,excluded=(['x_grid', 'y_grid', 'F','ghost_cells']))



def analytical(y,t):
  x, xdash, xddash = y
  dydt = [xdash, xddash, (2*np.pi*np.cos(2*np.pi*(t-x))*xddash /(np.sin(2*np.pi*(t-x))) ) -\
            (xdash * np.sin(2*np.pi*(t-x))**2 )+np.sin(2*np.pi*(t-x))\
         ]
  return dydt

initial_conditions_analytical = [initial_conditions_position_x[0], initial_conditions_velocity_x[0],initial_conditions_velocity_y[0] * np.sin(2*np.pi*(-initial_conditions_position_x[0]))]


Num_error = np.zeros(len(time), dtype = np.float)

""" Solving """

old = np.zeros(6 * no_of_particles, dtype=np.float)



old_analytical = np.zeros(6 * no_of_particles, dtype=np.float)


time_ana = np.arange(0, final_time, dt)
sol_analytical = odeint(analytical,initial_conditions_analytical,time_ana)
print(sol_analytical.shape)
""" Solver """

for time_index, t0 in enumerate(time):
    print("Computing for TimeIndex = ", time_index)
    t0 = time[time_index]
    if (time_index == time.size - 1):
        break
    t1 = time[time_index + 1]
    t = [t0, t1]
    if (time_index == 0):
        initial_conditions = initial_conditions
        #initial_conditions_analytical = initial_conditions_analytical
    else:
        initial_conditions = old
        #initial_conditions_analytical = old_analytical


    Ex_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_center,\
                                                  y_grid=y_top, F=Ex, ghost_cells = ghost_cells\
                                                )\
                          )

    Ey_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_right,\
                                                  y_grid=y_top, F=Ey, ghost_cells = ghost_cells\
                                                )\
                          )
    
    Ez_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_center,\
                                                  y_grid=y_top, F=Ez, ghost_cells = ghost_cells\
                                                )\
                          )

    Bx_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_right,\
                                                  y_grid=y_top, F=Bx, ghost_cells = ghost_cells\
                                                )\
                          )
    
    By_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_center,\
                                                  y_grid=y_top, F=By, ghost_cells = ghost_cells\
                                                )\
                          )

    Bz_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_right,\
                                                  y_grid=y_top, F=Bz, ghost_cells = ghost_cells\
                                                )\
                          )

    F_interpolated = np.concatenate([Ex_particle, Ey_particle, Ez_particle, Bx_particle, By_particle, Bz_particle], axis = 0)
    #print('Interpolated magnetic fields = ', F_interpolated[5,:])
    sol = mag_Verlet(initial_conditions, t, F_interpolated)
    
    

    #print('Analytical solution at current timestep = ',sol_analytical)
    Ex, Ey, Ez, Bx, By, Bz = fdtd(Ex, Ey, Ez, Bx, By, Bz, c, Lx, Ly, ghost_cells, Jx, Jy, Jz)
    
    
    #pl.contourf(x_right[ghost_cells:-ghost_cells], y_top[ghost_cells:-ghost_cells],Bz[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells], 100)
    #pl.show()
    #pl.clf()
    #print('Bz = ', Ey)
    #pl.contourf(x_right[ghost_cells:-ghost_cells],y_top[ghost_cells:-ghost_cells],Bz[ghost_cells:-ghost_cells,ghost_cells:-ghost_cells],100)
    #pl.plot(x_right[ghost_cells:-ghost_cells], Bz[5,ghost_cells:-ghost_cells], '-o')
    #pl.plot(x_center[ghost_cells:-ghost_cells],np.sin(2*np.pi*((time_index*dt)-x_center[ghost_cells:-ghost_cells])),'--')
    #pl.plot(x_right[ghost_cells:-ghost_cells], Ey[5,ghost_cells:-ghost_cells], '-o')
    #pl.plot(x_center[ghost_cells:-ghost_cells],np.sin(2*np.pi*((time_index*dt)-x_center[ghost_cells:-ghost_cells])),'--')
    #pl.ylim(-1,1)
    #pl.show()
    #pl.clf()
    
    for i in range(3 * no_of_particles):
        if (sol[i] >= right_boundary):
            sol[i] = sol[i] - Lx
        if (sol[i] <= left_boundary):
            sol[i] = sol[i] + Lx

    
    old = sol
    
    #old_analytical = sol_analytical[1, :]
    #sol_analytical = sol_analytical[1, :]
    #print(sol_analytical.shape)
    #print(sol_analytical)
    #print(sol_analytical[0])
    #Num_error[time_index] = abs ( old[0]-sol_analytical[0] )
    print('Numnerical x position = ', old[0])
    print('Analytical x position = ', sol_analytical[time_index + 1,0])
    zzz = input()
    #print('analytical x = ',sol_analytical)
    #h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'w')
    #h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=sol)
    #h5f.close()

    #h5f = h5py.File('analytical_all/solution_'+str(time_index)+'.h5', 'w')
    #h5f.create_dataset('analytical_all/solution_dataset_'+str(time_index), data=sol_analytical)
    #h5f.close()

    #print('shape is ', old_analytical.shape)

print(' the time scale is ', dt)
print('Error array  = ', Num_error)
