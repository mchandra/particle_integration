import numpy as np
import pylab as pl
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



def Interpolate(x, y, x_grid, y_grid, F):
    F_function = F.transpose()
    n      = len(F_function[0,:])-3           # number of zones

    x_zone = int(n * (x-x_grid[0]))             # indexing from zero itself
    y_zone = int(n * (y-y_grid[0]))

    b = np.matrix([ [0], [0], [0], [0] ])

    A = np.matrix( [ [1, x_grid[x_zone], y_grid[y_zone], x_grid[x_zone]*y_grid[y_zone] ], \
            [1, x_grid[x_zone], y_grid[y_zone + 1], x_grid[x_zone]*y_grid[y_zone + 1] ], \
            [1, x_grid[x_zone + 1], y_grid[y_zone], x_grid[x_zone + 1]*y_grid[y_zone] ], \
            [1, x_grid[x_zone + 1], y_grid[y_zone + 1], x_grid[x_zone + 1]*y_grid[y_zone + 1] ] ])



    point_to_calculated_for = np.matrix([ [1],[x], [y], [x*y] ])

    b = (la.inv(A)).transpose()*point_to_calculated_for


    Q11 = F_function[x_zone, y_zone]
    Q21 = F_function[x_zone + 1, y_zone]
    Q12 = F_function[x_zone, y_zone + 1]
    Q22 = F_function[x_zone + 1, y_zone + 1]

    Q = np.matrix([[Q11], [Q12], [Q21], [Q22] ])


    F_interpolated = b.transpose()*Q

    return F_interpolated









def error(a):
    nx = a  # number of zones not points
    ny = a  # number of zones not points

    dx = np.float(1 / (nx))
    dy = np.float(1 / (ny))

    x_center = np.linspace(-dx, 1 + dx, nx + 3, endpoint=True)
    y_center = np.linspace(-dy, 1 + dy, ny + 3, endpoint=True)

    Ez = np.zeros(((len(x_center)), (len(y_center))), dtype=np.float)
    Bx = np.zeros(((len(x_center)), (len(y_center) )), dtype=np.float)
    By = np.zeros(((len(x_center)), (len(y_center))), dtype=np.float)



    x_right = np.linspace(-dx / 2, 1 + 3*dx / 2, nx + 3, endpoint=True)
    y_top   = np.linspace(-dy / 2, 1 + 3*dy / 2, ny + 3, endpoint=True)

    # Conditions

    for i in range(nx + 1):
        for j in range(ny + 1):
            Ez[i + ghostcells][j + ghostcells] = np.sin(2 * np.pi * x_center[i+1] * y_center[j+1]) * np.cos(2 * np.pi * x_center[i+1] * y_center[j+1])
            Bx[i + ghostcells][j + ghostcells] = np.sin(2 * np.pi * x_center[i+1] * (y_top[j+1]-dy/2) ) * np.cos(2 * np.pi * x_center[i+1] * (y_top[j+1]-dy/2))
            By[i + ghostcells][j + ghostcells] = np.sin(2 * np.pi * (x_right[i+1]-dx/2) * y_center[j+1]) * np.cos(2 * np.pi * (x_right[i+1]-dx/2) * y_center[j+1])

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


    number_random_points = 30

    # Declaring random points

    x_random = np.random.rand(number_random_points)
    y_random = np.random.rand(number_random_points)



    Ez_at_random = np.zeros(number_random_points)
    Bx_at_random = np.zeros(number_random_points)
    By_at_random = np.zeros(number_random_points)


    # Calculating Interpolated values at the

    for i in range(number_random_points):
        Ez_at_random[i] = Interpolate( x_random[i] , y_random[i], x_center, y_center, Ez )
        Bx_at_random[i] = Interpolate(x_random[i] , y_random[i] , x_center, y_top, Bx)
        By_at_random[i] = Interpolate( x_random[i] , y_random[i], x_right, y_center, By )

    #Bx_at_random = BInterpolate(x_random , y_random , x_b, y_b, Bx[:-2, :-2])

    Ez_error = 0
    Bx_error = 0
    By_error = 0


    # error = sum( abs(Ez_at_random - np.sin(2*np.pi*x_random*y_random)*np.cos(2*np.pi*x_random*y_random)  ) ) /number_random_points
    Ez_error = sum(abs(Ez_at_random - np.sin(2 * np.pi * (x_random) * (y_random)) * np.cos( 2 * np.pi * (x_random) * (y_random)))) / number_random_points
    Bx_error = sum(abs(Bx_at_random - np.sin(2 * np.pi * (x_random) * (y_random - dy/2 )) * np.cos(2 * np.pi * (x_random) * (y_random - dy/2 ) )) ) / number_random_points
    By_error = sum(abs(By_at_random - np.sin(2 * np.pi * (x_random -dx/2) * (y_random)) * np.cos(2 * np.pi * (x_random-dx/2 ) * (y_random)))) / number_random_points

    return Ez_error,Bx_error, By_error


# Test Grid Sizes


N = np.arange(100,2000,100)
ErrorNEz = np.zeros(len(N), dtype=np.float)
ErrorNBx = np.zeros(len(N), dtype=np.float)
ErrorNBy = np.zeros(len(N), dtype=np.float)
for i in range(len(N)):
    ErrorNEz[i], ErrorNBx[i], ErrorNBy[i]  = error(N[i])
    print('Term = ', i)

# Plotting

pl.loglog(N, ErrorNEz, '-o', lw=3, label='$E_z$ ')
pl.legend()
pl.loglog(N, ErrorNBx, '-o', lw=3, label='$B_x$ ')
pl.legend()
pl.loglog(N, ErrorNBy, '-o', lw=5, label='$B_y$ ')
pl.legend()
pl.loglog(N, 150 * (N ** -1.999), '--', color='black', lw=2, label=' $O(N^{-2})$ ')
pl.legend()
pl.title('$\mathrm{Convergence\; plot}$ ')
pl.xlabel('$\mathrm{N}$')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
pl.show()
pl.clf()
