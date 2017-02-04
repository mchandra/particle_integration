from params import *
import numpy as np
from fields import *

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

    # J[y_current_zone, x_current_zone]+= (charge/(dx*dy))*velocity_required
    # return J
    return y_current_zone,x_current_zone,((charge/(dx*dy)))


charge_b0_depositor = np.vectorize(charge_b0_depositor, excluded=(['x_grid', 'y_grid', 'J','ghost_cells', 'Lx', 'Ly']))

# Example of usage Charge depositor

x = np.array([0.9, 0.1])
y = np.array([0.2, 0.6])

x_grid = np.array([-0.5, 0, 0.5, 1, 1.5])
y_grid = np.array([-0.5, 0, 0.5, 1, 1.5])

rho = np.matrix('0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0')

r = np.array(charge_b0_depositor(x = [x], y= [y], x_grid = x_grid, y_grid = y_grid, J = rho, ghost_cells = ghost_cells, Lx = Lx, Ly = Ly   ))

print(r)



"""Current Deposition for B0 splines (Have to vectorize)"""

def current_b0_depositor(x, y, velocity_required, x_grid, y_grid, J, ghost_cells, Lx, Ly):


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

    # J[y_current_zone, x_current_zone]+= (charge/(dx*dy))*velocity_required
    # return J
    return y_current_zone,x_current_zone,((charge/(dx*dy))*velocity_required)


current_b0_depositor = np.vectorize(current_b0_depositor, excluded=(['x_grid', 'y_grid', 'J','ghost_cells', 'Lx', 'Ly']))

# Example of usage Current depositor


x = np.array([0.2, 0.6])
y = np.array([0.2, 0.6])
v = np.array([5, 5])

x_grid = np.array([-0.5, 0, 0.5, 1, 1.5])
y_grid = np.array([-0.5, 0, 0.5, 1, 1.5])

J = np.matrix('0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0')

Jolo = np.array(current_b0_depositor(x = [x], y= [y], velocity_required = [v], x_grid = x_grid, y_grid = y_grid, J = J, ghost_cells = ghost_cells, Lx = Lx, Ly = Ly   ))

print(Jolo)


# J = current_b0_depositor(x = x, y= y, velocity_required = v, x_grid = x_grid, y_grid = y_grid,\
#                            J = J, ghost_cells = ghost_cells, Lx = Lx, Ly = Ly \
#                            )\
#
# print(J)



# def dcd(positions_n_plus_half ,velocites_n_plus_half, x_center_grid, y_center_grid,shape_function, no_of_particles, dx, dy):
#
#     x_right_grid = x_center_grid + dx
#     y_top_grid = y_center_grid + dy
#
#
#
#     Jx = shape_function(positions_n_plus_half[:no_of_particles], velocites_n_plus_half[:no_of_particles], x_right_grid, y_center_grid, dx, dy)
#     Jy = shape_function(positions_n_plus_half[no_of_particles:2*no_of_particles], velocites_n_plus_half[no_of_particles:2*no_of_particles], x_center_grid, y_top_grid, dx, dy)
#     Jz = shape_function(positions_n_plus_half[2*no_of_particles:3*no_of_particles], velocites_n_plus_half[2*no_of_particles:3*no_of_particles], x_center_grid, y_center_grid, dx, dy)
#
#     return Jx, Jy, Jz
