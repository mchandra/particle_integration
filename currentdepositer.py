from params import *
import numpy as np
from fields import *


# Ensure Particles have not crossed the domain boundaries in before current deposition

# Returns current at n+1/2 when x_n and x_n+1 are provided

# Ex = (x_right, y_center )
# Ey = (x_center, y_top   )
# Ez = (x_center, y_center)
# Bx = (x_center, y_top   )
# By = (x_right, y_center )
# Bz = (x_right, y_top    )
# rho = (x_center, y_top ) # Not needed here
# Jx = (x_right, y_center )
# Jy = (x_center, y_top )
# Jz = (x_center, y_center )

"""Charge Deposition for B0 splines (Have to vectorize)"""

def charge_b0_depositor(x, y, x_grid, y_grid, J, ghost_cells, Lx, Ly):


    n = (len(x_grid) - 1 - 2 * ghost_cells)  # number of zones

    dx = Lx/n
    dy = Ly/n

    x_zone = int(n * np.float(x - x_grid[0])/Lx)  # indexing from zero itself
    y_zone = int(n * np.float(y - y_grid[0])/Ly)

    if(abs(x-x_grid[x_zone])<abs(x-x_grid[x_zone + 1])):
        x_current_zone = x_zone
    else:
        x_current_zone = x_zone +1


    if(abs(y - y_grid[y_zone])<abs(y - y_grid[y_zone + 1])):
        y_current_zone = y_zone
    else:
        y_current_zone = y_zone +1

    return y_current_zone,x_current_zone,((charge/(dx*dy)))


charge_b0_depositor = np.vectorize(charge_b0_depositor, excluded=(['x_grid', 'y_grid', 'J','ghost_cells', 'Lx', 'Ly']))

# Example of usage Charge depositor

# x = np.array([0.9, 0.1])
# y = np.array([0.2, 0.6])
#
# x_grid = np.array([-0.5, 0, 0.5, 1, 1.5])
# y_grid = np.array([-0.5, 0, 0.5, 1, 1.5])
#
# rho = np.matrix('0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0')
#
# m,n,o = charge_b0_depositor(x = [x], y= [y], x_grid = x_grid, y_grid = y_grid, J = rho, ghost_cells = ghost_cells, Lx = Lx, Ly = Ly   )
#

# print(m)
# print(n)
# print(o)



"""Current Deposition for B0 splines (Have to vectorize)"""

def current_b0_depositor(charge, x, y, velocity_required, x_grid, y_grid, ghost_cells, Lx, Ly):


    n = (len(x_grid) - 1 - 2 * ghost_cells)  # number of zones

    dx = Lx/n
    dy = Ly/n
    x_zone = int(n * np.float(x - x_grid[0])/Lx)  # indexing from zero itself
    y_zone = int(n * np.float(y - y_grid[0])/Ly)

    if(abs(x-x_grid[x_zone])<abs(x-x_grid[x_zone + 1])):
        x_current_zone = x_zone
    else:
        x_current_zone = x_zone +1

    print('x_grid',len(x_grid))
    print('x', x)
    print('x_zone', x_zone)

    print('\n\n\n\n\n')
    #fasssssss

    print('y_grid',len(y_grid))
    print('y', y)
    print('y_zone', y_zone)




    if(abs(y - y_grid[y_zone])<abs(y - y_grid[y_zone + 1])):
        y_current_zone = y_zone
    else:
        y_current_zone = y_zone +1

    return y_current_zone,x_current_zone,((charge/(dx*dy))*velocity_required)


current_b0_depositor = np.vectorize(current_b0_depositor, excluded=(['charge','x_grid', 'y_grid', 'ghost_cells', 'Lx', 'Ly']))

# Example of usage Current depositor

# x = np.array([0.2, 0.6])
# y = np.array([0.2, 0.6])
# v = np.array([5, 5])
#
# x_grid = np.array([-0.5, 0, 0.5, 1, 1.5])
#
# y_grid = np.array([-0.5, 0, 0.5, 1, 1.5])
#
# J = np.matrix('0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0')
#
# i,j,k = np.array(current_b0_depositor(charge = 1, x = [x], y= [y], velocity_required = [v], x_grid = x_grid, y_grid = y_grid, ghost_cells = 1, Lx = 1, Ly = 1   ))
#
# # print(i,'\n \n')
#
# print('The i shape ', i.shape)
# print('The i  first element ', i[0,0])
#


# def fun(f):
#
#
#     x = np.array([0.2, 0.6])
#     y = np.array([0.2, 0.6])
#     v = np.array([5, 5])
#
#     x_grid = np.array([-0.5, 0, 0.5, 1, 1.5])
#
#     y_grid = np.array([-0.5, 0, 0.5, 1, 1.5])
#
#     e,b,c = f(charge = 1, x = [x], y= [y], velocity_required = [v], x_grid = x_grid, y_grid = y_grid, ghost_cells = 1, Lx = 1, Ly = 1   )
#     return e,b,c
#
# print(fun(current_b0_depositor))



def dcd(charge, no_of_particles, positions_plus_half ,velocities_plus_half, x_center_grid, y_center_grid,shape_function, ghost_cells, Lx, Ly, dx, dy):

    x_right_grid = x_center_grid + dx
    y_top_grid = y_center_grid + dy

    positions_x = positions_plus_half[:no_of_particles]
    positions_y = positions_plus_half[no_of_particles:2*no_of_particles]

    velocities_x = velocities_plus_half[:no_of_particles]
    velocities_y = velocities_plus_half[no_of_particles:2*no_of_particles]
    velocities_z = velocities_plus_half[2*no_of_particles:3*no_of_particles]


    Jx = np.zeros((len(x_center_grid), len(y_center_grid)), dtype=np.float)
    Jy = np.zeros((len(x_center_grid), len(y_center_grid)), dtype=np.float)
    Jz = np.zeros((len(x_center_grid), len(y_center_grid)), dtype=np.float)

    # print('fasfasfasfffff')



    # print('number of particle = ', no_of_particles)
    # print('length of positions given  ',len(positions_plus_half))
    # print('x',len(positions_x))
    # print('y',len(positions_y))
    # print('vx',len(velocities_x))
    print('y inside dcd function',positions_y)
    # print(shape_function(charge = charge,x = [positions_x], y = [positions_y], velocity_required = [velocities_x],x_grid = x_right_grid, y_grid = y_center_grid,\
    #                                                                         ghost_cells = ghost_cells,Lx =Lx, Ly= Ly))

    # print('Blashhhhhhhhhhhhh')

    Jx_x_indice, Jx_y_indices, Jx_values_at_these_indices = shape_function(charge = charge,x = [positions_x], y = [positions_y], velocity_required = [velocities_x],x_grid = x_right_grid, y_grid = y_center_grid,\
                                                                            ghost_cells = ghost_cells,Lx =Lx, Ly= Ly)


    for i in range(no_of_particles):
        Jx[Jx_x_indice[0,i], Jx_y_indices[0,i]] = Jx_values_at_these_indices[0,i]

    Jy_x_indice, Jy_y_indices, Jy_values_at_these_indices = shape_function(charge = charge,x = [positions_x], y = [positions_y], velocity_required = [velocities_y],x_grid = x_center_grid, y_grid = y_top_grid,\
                                                                            ghost_cells = ghost_cells,Lx =Lx, Ly= Ly)


    for i in range(no_of_particles):
        Jy[Jy_x_indice[0,i], Jy_y_indices[0,i]] = Jy_values_at_these_indices[0,i]



    Jz_x_indice, Jz_y_indices, Jz_values_at_these_indices = shape_function(charge = charge,x = [positions_x], y = [positions_y], velocity_required = [velocities_z],x_grid = x_center_grid, y_grid = y_center_grid,\
                                                                            ghost_cells = ghost_cells,Lx =Lx, Ly= Ly)


    for i in range(no_of_particles):
        Jy[Jz_x_indice[0,i], Jz_y_indices[0,i]] = Jz_values_at_these_indices[0,i]

    return Jx, Jy, Jz
