import numpy as np
import pylab as pl


"""  Plotting parameters  """
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

""" Plotting optimization"""



Nx = 100
Ny = 100

dx = float(1/Nx)
dy = float(1/Ny)



Ez = np.zeros((Nx,Ny),dtype = np.float)

Bx = np.zeros((Nx,Ny),dtype = np.float)

By = np.zeros((Nx,Ny),dtype = np.float)

x = np.zeros(Nx,dtype = np.float)
y = np.zeros(Ny,dtype = np.float)

x = np.linspace(0.5*dx,(Nx-0.5)*dx,Nx+1)
y = np.linspace(0.5*dx,(Ny-0.5)*dy,Ny+1)
spread = 8

c = 3.e8
dt = float(dx/(2*c))
dt_by_dx = dt/dx
dt_by_dy = dt/dy

max_iterations = 200






for time_index in range(max_iterations):
    print('Current time index is  ',time_index)
    
    for i in range(1,Nx):
        for j in range(1,Ny):
            
            Ez[i,j] = Ez[i,j] + ((dt_by_dx)*(By[i,j]-By[i-1,j])) - ((dt_by_dy)*(By[i,j]-By[i,j-1]))
            
    
            if(i==0 and j==0):
                Ez[i,j] = Ez[i,j] + ((dt_by_dx)*(By[i,j]-By[Nx+i-1,j])) - ((dt_by_dy)*(By[i,j]-By[i,Ny+j-1]))
            elif(i==0):
                Ez[i,j] = Ez[i,j] + ((dt_by_dx)*(By[i,j]-By[Nx+i-1,j])) - ((dt_by_dy)*(By[i,j]-By[i,j-1]))
            elif(j==0):
                Ez[i,j] = Ez[i,j] + ((dt_by_dx)*(By[i,j]-By[i-1,j])) - ((dt_by_dy)*(By[i,j]-By[i,Ny+j-1]))
    
    if(time_index==0):
        for i in range(Nx):
            for j in range(Ny):
                Ez[i,j] = np.exp(-((x[i]-0.5)**2+(y[j]-0.5)**2)/(2*spread**2))    
        
    for i in range(0,Nx-1):
        for j in range(0,Ny-1):
            
            Bx[i,j] = Bx[i,j] - ((dt_by_dy)*(Ez[i,j+1]-Ez[i,j]))
            By[i,j] = By[i,j] - ((dt_by_dx)*(Ez[i,j]-Ez[i+1,j])) 
                

            if(i==Nx-1 and j==Nx-1):
                Bx[i,j] = Bx[i,j] - ((dt_by_dy)*(Ez[i,j+1-Ny]-Ez[i,j]))
                By[i,j] = By[i,j] - ((dt_by_dx)*(Ez[i,j]-Ez[i+1-Nx,j]))                
            elif(i==Nx-1):
                Bx[i,j] = Bx[i,j] - ((dt_by_dy)*(Ez[i,j+1]-Ez[i,j]))
                By[i,j] = By[i,j] - ((dt_by_dx)*(Ez[i,j]-Ez[i+1-Nx,j]))                   
            elif(j==Ny-1):
                Bx[i,j] = Bx[i,j] - ((dt_by_dy)*(Ez[i,j+1-Ny]-Ez[i,j]))
                By[i,j] = By[i,j] - ((dt_by_dx)*(Ez[i,j]-Ez[i+1,j]))       

    #pl.contourf(Ez,100)
    #pl.title('Ez vs Time with Pulse')
    #pl.colorbar()
    #pl.savefig('images/point_mass' + '%04d'%time_index + '.png')

    #pl.clf()
    print(Ez)
