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
    dx = x_grid[x_zone + 1] - x_grid[x_zone]
    dy = y_grid[y_zone + 1] - y_grid[y_zone]

    E_interpolated = -(  (( (x_grid[x_zone+1] - x )*(y_grid[y_zone] - y ) )/( (dx) * (dy) ))*E[x_zone, y_zone + 1] +\
                     (( ( x - x_grid[x_zone] )*(y_grid[y_zone] - y ) )/( (dx) * (dy) )) * E[x_zone + 1, y_zone + 1] +\
                     (( (x_grid[x_zone+1] - x )*( y - y_grid[y_zone+1]) )/( (dx) * (dy) )) * E[x_zone, y_zone] +\
                     (( (x - x_grid[x_zone])*( y - y_grid[y_zone + 1]) )/( (dx) * (dy) )) * E[x_zone + 1, y_zone]  )


    return E_interpolated


# E_test = np.ones((3,3), dtype = np.float)
# E_test[:,1]*=2
# E_test[:,2]*=3
# E_test = np.matrix(E_test)
# print(E_test)
# print(E_test[1,0],E_test[1,1],E_test[1,2])
# print(E_test[0,1],E_test[1,1],E_test[2,1])
# print(    Interpolate( 0.5, 0.5, np.array([0,0.5,1]), np.array([0,0.5,1]) , E_test )    )


def error(a):


    nx = a # number of zones not points
    ny = a # number of zones not points

    x = np.linspace(0,1,nx+1,endpoint=True)
    x_plot = x[0:nx]


    y = np.linspace(0,1,ny+1, endpoint=True)
    y_plot = y[0:ny]


    Ez =  np.zeros(     (       (len(x)+2*ghostcells), ( len(y)+2*ghostcells )       ),dtype = np.float     )
    Bx =  np.zeros(     (       (len(x)+2*ghostcells), ( len(y)+2*ghostcells )       ),dtype = np.float     )
    By =  np.zeros(     (       (len(x)+2*ghostcells), ( len(y)+2*ghostcells )       ),dtype = np.float     )


    dx = np.float(1/(nx))
    dy = np.float(1/(ny))


    for i in range(nx+1):
        for j in range(ny+1):
            Ez[i + ghostcells][j + ghostcells] = np.exp(-x[i]**2-y[j]**2)#np.exp(-(y[j] - 0.5) ** 2 / (2 * spread ** 2))
            Bx[i + ghostcells][j + ghostcells] = np.exp(-((y[j] - 0.5) ** 2) / (2 * spread ** 2))
            By[i + ghostcells][j + ghostcells] = 0

#Ghost cells values
    Ez[0, :] = Ez[nx+1, :].copy()
    Ez[: , 0] = Ez[:, ny+1].copy()
    Ez[nx+1+ghostcells, :] = Ez[ghostcells, :].copy()
    Ez[ :, ny+1+ghostcells] = Ez[:, ghostcells].copy()

    Bx[0, :] = Bx[nx+1, :].copy()
    Bx[:, 0] = Bx[:, ny+1].copy()
    Bx[nx+1 + ghostcells, :] = Bx[ghostcells, :].copy()
    Bx[:, ny+1 + ghostcells] = Bx[:, ghostcells].copy()

    By[0, :] = By[nx+1, :].copy()
    By[:, 0] = By[:, ny+1].copy()
    By[nx+1 + ghostcells, :] = By[ghostcells, :].copy()
    By[:, ny+1 + ghostcells] = By[:, ghostcells].copy()

# Random points for error Testing

    number_random_points = 30

    x_random = np.random.rand(number_random_points)
    y_random = np.random.rand(number_random_points)


    Ez_at_random = np.random.rand(number_random_points)

    for i in range(number_random_points):
        Ez_at_random[i] = Interpolate( y_random[i], x_random[i], x, y,Ez[ghostcells:-ghostcells,ghostcells:-ghostcells].transpose())

    error = 0
    error = sum( abs(Ez_at_random - np.exp(-x_random**2-y_random**2)  ) ) /number_random_points
    return error


N = np.array([32, 64, 128, 256, 512,1024,2048,4096])
ErrorNEz = np.zeros(len(N), dtype=np.float)


for i in range(len(N)):
    ErrorNEz[i] = error(N[i])
    print('Term = ', i)


pl.loglog(N,ErrorNEz,'-o',lw =3,label = '$E_z$ ' )
pl.legend()
pl.loglog(N,1.5*(N**-1.999),'--',color = 'black',lw = 2,label = ' $O(N^{-2})$ ')
pl.legend()
pl.title('$\mathrm{Convergence\; plot}$ ')
pl.xlabel('$\mathrm{N}$')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')

pl.show()
pl.clf()


