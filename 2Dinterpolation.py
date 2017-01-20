import numpy as np
import pylab as pl
import numpy.linalg as la

import numpy.linalg as la

pl.rcParams['figure.figsize']   = 12, 7.5
pl.rcParams['lines.linewidth']  = 1.5
pl.rcParams['font.family']      = 'serif'
pl.rcParams['font.weight']      = 'bold'
pl.rcParams['font.size']        = 20
pl.rcParams['font.sans-serif']  = 'serif'
pl.rcParams['text.usetex']      = True
pl.rcParams['axes.linewidth']   = 1.5
pl.rcParams['axes.titlesize']   = 'medium'
pl.rcParams['axes.labelsize']   = 'medium'
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

# number of ghost cells

ghostcells = 1



def EInterpolate(x, y, x_grid, y_grid, E):
    E_function = E.transpose()
    n      = len(E_function[0,:])-1           # number of zones

    x_zone = int(n * x)             # indexing from zero itself
    y_zone = int(n * y)
    dx     = x_grid[x_zone + 1] - x_grid[x_zone]
    dy     = y_grid[y_zone + 1] - y_grid[y_zone]

    b = np.matrix([ [0], [0], [0], [0] ])

    A = np.matrix( [ [1, x_grid[x_zone], y_grid[y_zone], x_grid[x_zone]*y_grid[y_zone] ], \
            [1, x_grid[x_zone], y_grid[y_zone + 1], x_grid[x_zone]*y_grid[y_zone + 1] ], \
            [1, x_grid[x_zone + 1], y_grid[y_zone], x_grid[x_zone + 1]*y_grid[y_zone] ], \
            [1, x_grid[x_zone + 1], y_grid[y_zone + 1], x_grid[x_zone + 1]*y_grid[y_zone + 1] ] ])



    point_to_calculated_for = np.matrix([ [1],[x], [y], [x*y] ])

    b = (la.inv(A)).transpose()*point_to_calculated_for


    Q11 = E_function[x_zone, y_zone]
    Q21 = E_function[x_zone + 1, y_zone]
    Q12 = E_function[x_zone, y_zone + 1]
    Q22 = E_function[x_zone + 1, y_zone + 1]

    Q = np.matrix([[Q11], [Q12], [Q21], [Q22] ])


    E_interpolated = b.transpose()*Q

    return E_interpolated


# Correct

def BInterpolate(x, y, x_grid, y_grid, B):
    B_function = B.transpose()
    n      = len(B_function[0,:])-2         # number of zones

    x_zone = int(n * (x-x_grid[0]) )             # indexing from zero itself
    y_zone = int(n * (y-y_grid[0]))
    dx     = x_grid[x_zone + 1] - x_grid[x_zone]
    dy     = y_grid[y_zone + 1] - y_grid[y_zone]

    #print('Magnetic Field = ', B)
    b = np.matrix([ [0], [0], [0], [0] ])

    A = np.matrix( [ [1, x_grid[x_zone], y_grid[y_zone], x_grid[x_zone]*y_grid[y_zone] ], \
            [1, x_grid[x_zone], y_grid[y_zone + 1], x_grid[x_zone]*y_grid[y_zone + 1] ], \
            [1, x_grid[x_zone + 1], y_grid[y_zone], x_grid[x_zone + 1]*y_grid[y_zone] ], \
            [1, x_grid[x_zone + 1], y_grid[y_zone + 1], x_grid[x_zone + 1]*y_grid[y_zone + 1] ] ])



    point_to_calculated_for = np.matrix([ [1],[x], [y], [x*y] ])

    b = (la.inv(A)).transpose()*point_to_calculated_for


    Q11 = B_function[x_zone, y_zone]
    Q21 = B_function[x_zone + 1, y_zone]
    Q12 = B_function[x_zone, y_zone + 1]
    Q22 = B_function[x_zone + 1, y_zone + 1]

    Q = np.matrix([[Q11], [Q12], [Q21], [Q22] ])


    B_interpolated = b.transpose()*Q

    # print('Inside Function n is  = ', n)
    # print('x is =', x)
    # print('y is =', y)
    # print('Input x grid  = ', x_grid)
    # print('Input y grid  = ', y_grid)
    # print('x_zone  = ', x_zone)
    # print('y_zone  = ', y_zone)
    # print('Input Magnetic Field  = ', B)
    #
    # print('Inside Function Interpolated Magenetic Field = ', B_interpolated)
    return B_interpolated





#
# B_test = np.ones((4,4), dtype = np.float)
# B_test[:,1]*=2
# B_test[:,2]*=3
# print(B_test)
# print(B_test[1,0],B_test[1,1],B_test[1,2],B_test[1])
# print(B_test[0,1],B_test[1,1],B_test[2,1])
# print(    BInterpolate( 0.5, 0.25, np.array([-0.25, 0.25, 0.75, 1.25]), np.array([-0.25, 0.25, 0.75, 1.25]) , B_test )    )
#
#





# E_test = np.ones((3,3), dtype = np.float)
# E_test[:,1]*=2
# E_test[:,2]*=3
# print(E_test)
# print(E_test[1,0],E_test[1,1],E_test[1,2])
# print(E_test[0,1],E_test[1,1],E_test[2,1])
# print(    EInterpolate( 0.25, 0.75, np.array([0,0.5,1]), np.array([0,0.5,1]) , E_test )    )



def error(a):
    nx = a  # number of zones not points
    ny = a  # number of zones not points

    x = np.linspace(0, 1, nx + 1, endpoint=True)
    y = np.linspace(0, 1, ny + 1, endpoint=True)

    Ez = np.zeros(((len(x) + 2 * ghostcells), (len(y) + 2 * ghostcells)), dtype=np.float)
    Bx = np.zeros(((len(x) + 2 * ghostcells), (len(y) + 2 * ghostcells)), dtype=np.float)
    By = np.zeros(((len(x) + 2 * ghostcells), (len(y) + 2 * ghostcells)), dtype=np.float)

    dx = np.float(1 / (nx))
    dy = np.float(1 / (ny))

    x_b = np.linspace(-dx / 2, 1 + dx / 2, nx + 2, endpoint=True)
    y_b = np.linspace(-dx / 2, 1 + dx / 2, ny + 2, endpoint=True)

    # Conditions

    for i in range(nx + 1):
        for j in range(ny + 1):
            Ez[i + ghostcells][j + ghostcells] = np.sin(2 * np.pi * x[i] * y[j]) * np.cos(2 * np.pi * x[i] * y[j])
            Bx[i + ghostcells][j + ghostcells] = np.sin(2 * np.pi * x[i] * y[j]) * np.cos(2 * np.pi * x[i] * y[j])
            By[i + ghostcells][j + ghostcells] = np.sin(2 * np.pi * (x_b[i + 1]) * y_b[j + 1]) * np.cos(
                2 * np.pi * x_b[i + 1] * y_b[j + 1])

            # Ghost cells values copying

    Ez[0, :] = Ez[nx + 1, :].copy()
    Ez[:, 0] = Ez[:, ny + 1].copy()
    Ez[nx + 1 + ghostcells, :] = Ez[ghostcells, :].copy()
    Ez[:, ny + 1 + ghostcells] = Ez[:, ghostcells].copy()

    Bx[0, :] = Bx[nx + 1, :].copy()
    Bx[:, 0] = Bx[:, ny + 1].copy()
    Bx[nx + 1 + ghostcells, :] = Bx[ghostcells, :].copy()
    Bx[:, ny + 1 + ghostcells] = Bx[:, ghostcells].copy()

    By[0, :] = By[nx + 1, :].copy()
    By[:, 0] = By[:, ny + 1].copy()
    By[nx + 1 + ghostcells, :] = By[ghostcells, :].copy()
    By[:, ny + 1 + ghostcells] = By[:, ghostcells].copy()

    # Random points for error Testing

    number_random_points = 50

    # Declaring random points

    x_random = np.random.rand(number_random_points)
    y_random = np.random.rand(number_random_points)

    # x_random = np.array(0.1)
    # y_random = np.array(0.1)
    x_random_b = x_random
    y_random_b = y_random

    # x_random = np.array([0.50,0.50,0.50])
    # y_random = np.array([0.25,0.5,0.75])
    Ez_at_random = np.zeros(number_random_points)
    Bx_at_random = np.zeros(number_random_points)
    By_at_random = np.zeros(number_random_points)

    #print('Physical grid = ', x)
    # Calculating Interpolated values at the

    for i in range(number_random_points):
        # Ez_at_random[i] = EInterpolate( y_random[i], x_random[i], x, y, Ez[ghostcells:-ghostcells,ghostcells:-ghostcells].transpose() )
        Bx_at_random[i] = BInterpolate(x_random[i] , y_random[i] , x_b, y_b, Bx[:-1, :-1])
        # By_at_random[i] = BInterpolate( y_random[i] - dy/2, x_random[i] - dx/2, x_b, y_b, By[:-ghostcells, ghostcells:-ghostcells].transpose() )

    #Bx_at_random = BInterpolate(x_random , y_random , x_b, y_b, Bx[:-2, :-2])
    error = 0
    # error = sum( abs(Ez_at_random - np.sin(2*np.pi*x_random*y_random)*np.cos(2*np.pi*x_random*y_random)  ) ) /number_random_points
    error = sum(abs(Bx_at_random - np.sin(2 * np.pi * (x_random-dx/2) * (y_random-dy/2)) * np.cos(2 * np.pi * (x_random-dx/2) * (y_random-dy/2)) )) / number_random_points
    #print('expected Magnetic field = ',(x_random-dx/2)**2+(y_random-dy/2)**2 )
    #print('dx = ', dx)
    #print('error = ',error)
    return error


# Test Grid Sizes

#N = np.array([32,64,128,256,512])
N = np.arange(30,600,30)
ErrorNEz = np.zeros(len(N), dtype=np.float)

for i in range(len(N)):
    ErrorNEz[i] = error(N[i])
    print('Term = ', i)

# Plotting

pl.loglog(N, ErrorNEz, '-o', lw=3, label='$E_z$ ')
pl.legend()
pl.loglog(N, 150 * (N ** -1.999), '--', color='black', lw=2, label=' $O(N^{-2})$ ')
pl.legend()
pl.title('$\mathrm{Convergence\; plot}$ ')
pl.xlabel('$\mathrm{N}$')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
pl.show()
pl.clf()
