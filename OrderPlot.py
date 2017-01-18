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
    return np.float(nodal_point/(local_domain_nodes-1))

def shape(localXi,local_array): 
    #print('The Xi for the term is ')
    
    local_domain_nodes = len(local_array)
    localXi_matrix = [[ 0 ] for i in range(local_domain_nodes)]
    u = np.matrix(local_array)
    
    for i in range(local_domain_nodes):
        localXi_matrix[i] = localXi**i
    #print('local Xi vector matrix is ',localXi_matrix)
    c = np.matrix( [[ 0 ] for i in range(local_domain_nodes)] )
    A = np.zeros( (local_domain_nodes , local_domain_nodes ) ,dtype = np.float)
    
    for i in range(local_domain_nodes):
        for j in range(local_domain_nodes):
            A[i,j] = getXi(i, local_domain_nodes)**j
    #print('A is  ',A)
    c = np.matrix(la.inv(A))*u.transpose()
    #print('c is  is ',c)
    #print('The answer returned is  ', np.matrix(localXi_matrix)*np.matrix(c))
    return np.matrix(localXi_matrix)*np.matrix(c)


def error(Nx,order):
    number_of_ghost_points = 0
    x1 = np.linspace(0, 1, Nx+1)
    x2 = np.linspace(0, 1, Nx+1)

    Electric_field = np.sin(4*np.pi*x1)

    number_of_random_points = 10

    random_points = np.random.rand(number_of_random_points)
    Xi_random_points = np.zeros(number_of_random_points, dtype = np.float)
    zones_random_points = np.zeros(number_of_random_points, dtype = np.float)


    actual_electric_fields = np.sin(4*np.pi*random_points)

    interpolated_electric_fields = np.zeros(len(random_points), dtype = np.float)

    order_of_interpolation = order


    Xi_random_points =   ( (  (random_points*Nx)%(order_of_interpolation)  )/order_of_interpolation  )
    for i in range(number_of_random_points):
        zones_random_points[i] =     np.int(   (   (random_points[i]*Nx)/order_of_interpolation   )  )+1


    if(np.int( (zones_random_points[i]-1)*(order_of_interpolation) )  + order_of_interpolation+1 > Nx):
        number_of_ghost_points = np.int( (zones_random_points[i])*(order_of_interpolation) )  +1 - Nx

        Electric_field = np.concatenate([Electric_field,Electric_field[0:number_of_ghost_points]],axis = 0)


    for i in range(number_of_random_points):
        interpolated_electric_fields[i] = shape(Xi_random_points[i],Electric_field[  np.int( (zones_random_points[i]-1)*(order_of_interpolation) )   :     np.int( (zones_random_points[i]-1)*(order_of_interpolation) )  + order_of_interpolation+1      ] )



    x_refined = np.linspace(0,1,101)
    E_refined = np.sin(4*np.pi*x_refined)

    absError = sum(abs(interpolated_electric_fields-actual_electric_fields))/number_of_random_points
    
    return absError




#print(shape(0.35,np.array([1,2])))



N = np.array( [ 32,64,128,256,512,1024, 2048, 4096, 8192  ] )
ErrorNEOne = np.zeros((  len(N) )  ,dtype = np.float)
ErrorNEThree = np.zeros((  len(N) )  ,dtype = np.float)
ErrorNEFive = np.zeros((  len(N) )  ,dtype = np.float)
ErrorNESeven = np.zeros((  len(N) )  ,dtype = np.float)
ErrorNENine = np.zeros((  len(N) )  ,dtype = np.float)
for i in range(len(N)):

    ErrorNEOne[i] = error(N[i],1)
    
for i in range(len(N)):

    ErrorNEThree[i] = error(N[i],3)
    
for i in range(len(N)):

    ErrorNEFive[i] = error(N[i],5)
    
for i in range(len(N)):

    ErrorNESeven[i] = error(N[i],7)
    
for i in range(len(N)):

    ErrorNENine[i] = error(N[i],9)
    

pl.loglog(N,ErrorNEOne,'-o',label = '$Order= '+str(1)+'$' )
pl.legend().draggable()
pl.loglog(N,ErrorNEThree,'-o',label = '$Order= '+str(3)+'$' )
pl.legend().draggable()
pl.loglog(N,ErrorNEFive,'-o',label = '$Order= '+str(5)+'$' )
pl.legend().draggable()
pl.loglog(N,ErrorNESeven,'-o',label = '$Order= '+str(7)+'$' )
pl.legend().draggable()
pl.loglog(N,ErrorNENine,'-o',label = '$Order= '+str(9)+'$' )
pl.legend().draggable()
#pl.loglog(N,ErrorNE[:,1],'-o',lw =1,label = '$Order= '+str(3)+'$' )
#pl.legend()
#pl.loglog(N,ErrorNE[:,2],lw =1,label = '$Order= '+str(5)+'$' )
#pl.legend()
#pl.loglog(N,ErrorNE[:,3],'-o',lw =1,label = '$Order= '+str(7)+'$' )
#pl.legend()
#pl.loglog(N,ErrorNE[:,4],lw =1,label = '$Order= '+str(9)+'$' )
#pl.legend()
pl.loglog(N,0.15*(N**-.999),'--',color = 'black',label = ' $O(N^{-1})$ ')
pl.legend().draggable()
pl.title('$\mathrm{Convergence\; plot}$ ')
pl.xlabel('$\mathrm{N}$')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
pl.show()
pl.clf()
