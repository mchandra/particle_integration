from params import *
from scipy.integrate import odeint


#Equations are:
# dBz/dt = - dEy/dx
# dEy/dt = - dBz/dx

# x, E in time n, n+1.............
# v, B in time n +0.5, n + 1.5 .....
# for Boris Algorithm : x(n+1) = x(n) + v(n+0.5)dt
#  v(n+1.5) = v(n + 0.5) + fields(E(n+1), B(avg(n+1.5,n+0.5)))


# For analytical comparision
# for x start from x(n+1) and for vx start from time (n+1)

""" Setting number of particles and other parameters"""

no_of_particles = 50

Nx = 100
Ny = 100

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


""" Initial conditions for positions """
# At n = 0
left_boundary = 0
right_boundary = Lx
length_of_box_x = right_boundary - left_boundary
initial_conditions_position_x = left_boundary + length_of_box_x * np.random.rand(no_of_particles)
# At n = 0
bottom_boundary = 0
top_boundary = Ly
length_of_box_y = length_of_box_x
initial_conditions_position_y = bottom_boundary + length_of_box_y * np.random.rand(no_of_particles)
# At n = 0
back_boundary = 0
front_boundary = Ly
length_of_box_z = length_of_box_x
initial_conditions_position_z = back_boundary + length_of_box_z * np.random.rand(no_of_particles)

""" Setting velocities according to maxwellian distribution """
# At n = 0.5
initial_conditions_velocity_x = np.random.rand(no_of_particles)*0.2
initial_conditions_velocity_y = np.random.rand(no_of_particles)*0.2
initial_conditions_velocity_z = np.zeros(no_of_particles, dtype=np.float)

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

# mass and charge of particle:

me = 1
charge = 1

""" Writing the spatial grids as a two dimension matrix for vectorization purposes """

X_center_physical, Y_center_physical = np.meshgrid( x_center[ghost_cells:-ghost_cells],\
                                                    y_center[ghost_cells:-ghost_cells]\
                                                  )

""" Discretizing time and making sure scaling is done right """

# box_crossing_time_scale = length_of_box_x / np.max(initial_conditions_velocity_x)

final_time = 5
dt = np.float(dx / (2 * c))
time = np.arange(0, final_time, dt)

""" Writing the offset spatial grids and indices as a two dimension matrix for vectorization purposes """

X_right_physical, Y_top_physical = np.meshgrid( x_right[ghost_cells:-ghost_cells],\
                                                y_top[ghost_cells:-ghost_cells]\
                                              )

I, J = np.meshgrid( range(ghost_cells, len(x_center) - ghost_cells),\
                    range(ghost_cells, len(y_center) - ghost_cells)\
                  )
# Initializing the non relevant fields:

Ey[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = np.sin(2*np.pi*(-X_right_physical))
Bz[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = np.sin(2*np.pi*((dt/2)-X_right_physical))

#Bz[ghost_cells:-ghost_cells, ghost_cells:-ghost_cells] = 20


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

def analytical(Y,t):
  x, y, vx, vy = Y
  dydt = [ vx, vy, vy*np.sin(2*np.pi*(t + dt - x)), (1-vx)*np.sin(2*np.pi*(t + dt - x)) ]
  return dydt

position_numerical = np.zeros((len(time),2), dtype = np.float)
velocity_numerical = np.zeros((len(time),2), dtype = np.float)

Num_error = np.zeros((len(time),2), dtype = np.float)


""" Solving """

old = np.zeros(6 * no_of_particles, dtype=np.float)

""" Solver """

for time_index, t0 in enumerate(time):
    print("Computing for TimeIndex = ", time_index)
    # print('\n \n')
    t0 = time[time_index]
    if (time_index == time.size - 1):
        break
    t1 = time[time_index + 1]
    t = [t0, t1]
    if (time_index == 0):
        initial_conditions = initial_conditions

    else:
        initial_conditions = old


    Ex_updated, Ey_updated, Ez_updated, Bx_updated, By_updated, Bz_updated = fdtd(Ex, Ey, Ez, Bx, By, Bz, c, Lx, Ly, ghost_cells, Jx, Jy, Jz)


    Ex_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_center,\
                                                  y_grid=y_top, F=Ex_updated, ghost_cells = ghost_cells, Lx = Lx, Ly = Ly\
                                                )\
                          )

    Ey_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_right,\
                                                  y_grid=y_top, F=Ey_updated, ghost_cells = ghost_cells, Lx = Lx, Ly = Ly\
                                                )\
                          )

    Ez_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_center,\
                                                  y_grid=y_top, F=Ez_updated, ghost_cells = ghost_cells, Lx = Lx, Ly = Ly\
                                                )\
                          )

    Bx_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_right,\
                                                  y_grid=y_top, F=((Bx+Bx_updated)/2), ghost_cells = ghost_cells, Lx = Lx, Ly = Ly\
                                                )\
                          )

    By_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_center,\
                                                  y_grid=y_top, F=((By+By_updated)/2), ghost_cells = ghost_cells, Lx = Lx, Ly = Ly\
                                                )\
                          )

    Bz_particle = np.array( bilinear_interpolate( x=[initial_conditions[:no_of_particles]], y=[initial_conditions[no_of_particles:2*no_of_particles]], x_grid=x_right,\
                                                  y_grid=y_top, F=((Bz+Bz_updated)/2), ghost_cells = ghost_cells, Lx = Lx, Ly = Ly\
                                                )\
                          )

    F_interpolated = np.concatenate([Ex_particle, Ey_particle, Ez_particle, Bx_particle, By_particle, Bz_particle], axis = 0)

    sol = mag_Verlet(initial_conditions, t, F_interpolated)

    indices_right = np.where(sol[0:3*no_of_particles]>right_boundary)
    indices_left = np.where(sol[0:3*no_of_particles]<left_boundary)
    sol[indices_right] -= Lx
    sol[indices_left] += Lx

    Ex, Ey, Ez, Bx, By, Bz= Ex_updated, Ey_updated, Ez_updated, Bx_updated, By_updated, Bz_updated

    old = sol

    h5f = h5py.File('solution_all/solution_'+str(time_index)+'.h5', 'w')
    h5f.create_dataset('solution_all/solution_dataset_'+str(time_index), data=sol)
    h5f.close()
