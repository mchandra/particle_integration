import numpy as np
import pylab as pl
import numpy.linalg as la

pl.rcParams['figure.figsize'] = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family'] = 'serif'
pl.rcParams['font.weight'] = 'bold'
pl.rcParams['font.size'] = 20
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex'] = True
pl.rcParams['axes.linewidth'] = 1.5
pl.rcParams['axes.titlesize'] = 'medium'
pl.rcParams['axes.labelsize'] = 'medium'

pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.major.pad'] = 8
pl.rcParams['xtick.minor.pad'] = 8
pl.rcParams['xtick.color'] = 'k'
pl.rcParams['xtick.labelsize'] = 'medium'
pl.rcParams['xtick.direction'] = 'in'

pl.rcParams['ytick.major.size'] = 8
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.major.pad'] = 8
pl.rcParams['ytick.minor.pad'] = 8
pl.rcParams['ytick.color'] = 'k'
pl.rcParams['ytick.labelsize'] = 'medium'
pl.rcParams['ytick.direction'] = 'in'



def sumsum(a):
    return sum(sum(abs(a) ) )


spread = 0.1
ghostcells = 1


def Interpolate(y, x, x_grid, y_grid, E):
    n = len(E[0,:])-2
    x_zone = int(n * x)             # indexing from zero itself
    y_zone = int(n * y)
    print('n = ', n)
    print('printing x = ',x)
    print('printing y = ', y)
    print('printing x_grid',x_grid)
    print('printing y_grid', y_grid)
    print('printing x zone',x_zone)
    print('printing y zone', y_zone)
    print('Electric field  = ',E)
    dx = x_grid[x_zone + 1] - x_grid[x_zone]
    dy = y_grid[y_zone + 1] - y_grid[y_zone]

    E_interpolated = -(  (( (x_grid[x_zone+1] - x )*(y_grid[y_zone] - y ) )/( (dx) * (dy) ))*E[x_zone, y_zone + 1] +\
                     (( ( x - x_grid[x_zone] )*(y_grid[y_zone] - y ) )/( (dx) * (dy) )) * E[x_zone + 1, y_zone + 1] +\
                     (( (x_grid[x_zone+1] - x )*( y - y_grid[y_zone+1]) )/( (dx) * (dy) )) * E[x_zone, y_zone] +\
                     (( (x - x_grid[x_zone])*( y - y_grid[y_zone + 1]) )/( (dx) * (dy) )) * E[x_zone + 1, y_zone]  )


    print('Interpolated answer = ', E_interpolated)

    return E_interpolated


# E_test = np.ones((3,3), dtype = np.float)
# E_test[:,1]*=2
# E_test[:,2]*=3
# E_test = np.matrix(E_test)
# print(E_test)
# print(E_test[1,0],E_test[1,1],E_test[1,2])
# print(E_test[0,1],E_test[1,1],E_test[2,1])
# print(    Interpolate( 0.5, 0.5, np.array([0,0.5,1]), np.array([0,0.5,1]) , E_test )    )
