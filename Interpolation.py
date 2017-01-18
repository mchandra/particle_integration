import numpy as np
import pylab as pl
import numpy.linalg as la


pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
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




def getXi(nodal_point, local_domain_nodes):
    return np.float(nodal_point/local_domain_nodes)



def shape(localXi,local_array): 
    
    local_domain_nodes = len(local_array)
    localXi_matrix = [[ 0 ] for i in range(local_domain_nodes)]
    u = np.matrix(local_array)
    
    for i in range(local_domain_nodes):
        localXi_matrix[i] = localXi**i
        
    c = np.matrix( [[ 0 ] for i in range(local_domain_nodes)] )
    A = np.zeros( (local_domain_nodes , local_domain_nodes ) ,dtype = np.float)
    
    for i in range(local_domain_nodes):
        for j in range(local_domain_nodes):
            A[i,j] = getXi(i, local_domain_nodes)**j
    
    c = np.matrix(la.inv(A))*u.transpose()

    return np.matrix(localXi_matrix)*np.matrix(c)


Nx = 1000

x1 = np.linspace(0, 1, Nx+1)
x2 = np.linspace(0, 1, Nx+1)

Electric_field = np.sin(4*np.pi*x1)
Magnetic_field = np.sin(4*np.pi*x2)


number_of_random_points = 20

random_points = np.random.rand(number_of_random_points)
Xi_random_points = np.zeros(number_of_random_points, dtype = np.float)
zones_random_points = np.zeros(number_of_random_points, dtype = np.float)


actual_electric_fields = np.sin(4*np.pi*random_points)

interpolated_electric_fields = np.zeros(len(random_points), dtype = np.float)

order_of_interpolation = 2



Xi_random_points =   ( (  (random_points*Nx)%(order_of_interpolation)  )/order_of_interpolation  )
for i in range(number_of_random_points):
    zones_random_points[i] =     np.int(   (   (random_points[i]*Nx)/order_of_interpolation   )  )



for i in range(number_of_random_points):
    interpolated_electric_fields[i] = shape(Xi_random_points[i],Electric_field[  np.int( (zones_random_points[i]-1)*(order_of_interpolation+1) -1)   :     np.int( (zones_random_points[i]-1)*(order_of_interpolation+1) -1)  + order_of_interpolation+1      ] )


x_refined = np.linspace(0,1,101)
E_refined = np.sin(4*np.pi*x_refined)

pl.plot(random_points, interpolated_electric_fields,'*')
pl.plot(random_points, actual_electric_fields,'o')
pl.plot(x_refined,E_refined)
pl.show()
pl.clf()


