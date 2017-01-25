import numpy as np
import pylab as pl

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
pl.rcParams['text.latex.unicode'] = True
pl.rcParams['xtick.major.size'] = 8     
pl.rcParams['xtick.minor.size'] = 4     
pl.rcParams['xtick.major.pad']  = 8     
pl.rcParams['xtick.minor.pad']  = 8     
pl.rcParams['xtick.color']      = 'k'     
pl.rcParams['xtick.labelsize']  = 'medium'
pl.rcParams['xtick.direction']  = 'in'    
pl.rcParams['image.aspect']  = 'equal'    
pl.rcParams['ytick.major.size'] = 8     
pl.rcParams['ytick.minor.size'] = 4     
pl.rcParams['ytick.major.pad']  = 8     
pl.rcParams['ytick.minor.pad']  = 8     
pl.rcParams['ytick.color']      = 'k'     
pl.rcParams['ytick.labelsize']  = 'medium'
pl.rcParams['ytick.direction']  = 'in'    



# dEz/dt = dBy/dx - dBx/dy
# dBx/dt = -dEz/dy
# dBy/dt = +dEz/dx
# div_B = dBx/dx + dBy/dy


spread = 0.1


def gauss(zvar,time):
    
    return np.exp(  -  (zvar-0.5 )**2/(2*spread**2))  #(0.5+time*dt) - int(0.5+time*dt)

ghostcells = 1
Nx = 100
Ny = 100


absSumErrorEz = 0
absSumErrorBx = 0
absSumErrorBy = 0



x = np.linspace(0,1,Nx+1)
x_plot = x[0:Nx]


y = np.linspace(0,1,Ny+1)
y_plot = y[0:Ny]


Ez =  np.zeros(     (       (Nx+2*ghostcells), ( Ny+2*ghostcells )       ),dtype = np.float     )
Bx =  np.zeros(     (       (Nx+2*ghostcells), ( Ny+2*ghostcells )       ),dtype = np.float     )
By =  np.zeros(     (       (Nx+2*ghostcells), ( Ny+2*ghostcells )       ),dtype = np.float     )


div_B = np.zeros(  ( (Nx +2*ghostcells ), ( Ny +2*ghostcells )   ),dtype = np.float)
c = 1
div_B_ini = np.zeros(  ( (Nx +2*ghostcells ), ( Ny +2*ghostcells )   ),dtype = np.float)

dx = np.float(1/(Nx))
dy = np.float(1/(Ny))
dt = np.float(dx/(2*c))


max_iterations = np.int(4/(dt))

dt_by_dx = dt/(dx)
dt_by_dy = dt/(dy)

for time_index in range(max_iterations):
    print('Time Index = ', time_index)
    if(time_index==0):
        for i in range(Nx):                
            for j in range(Ny):
                Ez[i+ghostcells, j+ghostcells] = np.exp(- (     (x[i]-0.5)**2   +     (y[j]-0.5)**2      )     /(2*spread**2))
                Bx[ i+ghostcells, j+ghostcells] = np.exp(-(y[j]-0.5)**2/(2*spread**2))
                By[i+ghostcells, j+ghostcells] = np.exp(-(x[i]-0.5)**2/(2*spread**2))
                
        Bx[0, :] = Bx[Nx, :]
        Bx[:, 0] = Bx[:, Ny]
        Bx[Nx+ghostcells, :] = Bx[ghostcells,:]
        Bx[:, Ny+ghostcells] = Bx[:, ghostcells]
        
        By[0, :] = By[Nx, :]
        By[:, 0] = By[:, Ny]
        By[Nx+ghostcells, :] = By[ghostcells,:]
        By[:, Ny+ghostcells] = By[:, ghostcells]                
        for i in range(ghostcells,Nx+ghostcells):
            for j in range(ghostcells, Ny + ghostcells):
            
                div_B_ini[i, j] = (Bx[i, j]-Bx[i-1,j])/(dx) +  (By[i, j]-By[i, j- 1])/(dy)        
       
       
        pl.figure(figsize = (13,10) )
        pl.contourf(x_plot, y_plot ,div_B_ini[1:Nx+1,1:Ny+1],100)
        pl.title( r'$\nabla \cdot \mathbf{B}$')
        pl.xlabel('$x$ ')
        pl.ylabel('$y$')
        pl.colorbar()
        #pl.grid(True)
        #pl.ylim([0,1])
        #pl.xlim([0,1])
        #pl.axes().set_aspect('equal', 'datalim')
        #pl.draw()
        pl.savefig('div/point_mass' + '%04d'%0 + '.png')
        #pl.autoscale()
        #pl.show()
        pl.clf()        
        pl.close()
        
    for i in range(ghostcells, Nx + ghostcells):
        for j in range(ghostcells, Ny + ghostcells):
            
            Ez[i, j] = Ez[i, j] + (dt_by_dx*(By[i, j]-By[i-1, j])) -  (dt_by_dy*(Bx[i, j]-Bx[i,  j - 1 ]))


    
    Ez[0, :] = Ez[Nx, :]
    Ez[: , 0] = Ez[:, Ny]
    
    Ez[Nx+ghostcells, :] = Ez[ghostcells, :]
    Ez[ :, Ny+ghostcells] = Ez[:, ghostcells]
    
    
    
    
    for i in range(ghostcells,Nx+ghostcells):
        for j in range(ghostcells, Ny + ghostcells):
            
            Bx[i, j] = Bx[i, j] - (dt_by_dy*(Ez[i, j+1] - Ez[i, j]))
            
            By[i, j] = By[i, j] + (dt_by_dx*(Ez[i+1, j] - Ez[i, j]))
    
    
    
    Bx[0, :] = Bx[Nx, :]
    Bx[:, 0] = Bx[:, Ny]
    Bx[Nx+ghostcells, :] = Bx[ghostcells,:]
    Bx[:, Ny+ghostcells] = Bx[:, ghostcells]
    
    By[0, :] = By[Nx, :]
    By[:, 0] = By[:, Ny]
    By[Nx+ghostcells, :] = By[ghostcells,:]
    By[:, Ny+ghostcells] = By[:, ghostcells]
    
    for i in range(ghostcells,Nx+ghostcells):
        for j in range(ghostcells, Ny + ghostcells):
            #print((Bx[i, j]-Bx[i-1,j]))
           # print((By[i, j]-By[i, j- 1]))
            div_B[i, j] = (Bx[i+1, j]-Bx[i,j])/(dx) +  (By[i, j+1]-By[i, j])/(dy)

    div_B[0, :] = div_B[Nx, :]
    div_B[:, 0] = div_B[:, Ny]
    div_B[Nx+ghostcells, :] = div_B[ghostcells,:]
    div_B[:, Ny+ghostcells] = div_B[:, ghostcells]
    
    tim = time_index+1
    
    pl.figure(figsize = (13,10) )
    pl.contourf(x_plot, y_plot ,Ez[1:Nx+1,1:Ny+1] ,100)
    pl.title(' $E_z(x, y)$')
    pl.xlabel('$x$ ')
    pl.ylabel('$y$')
    pl.colorbar()
    #pl.ylim([0,1])
    #pl.xlim([0,1])    
    #pl.grid(True)
    #pl.axes().set_aspect('equal', 'datalim')
    pl.savefig('images/point_mass' + '%04d'%time_index + '.png')
    pl.clf()            
    pl.close()
    
    
    pl.figure(figsize = (13,10) )
    pl.contourf(x_plot, y_plot ,By[1:Nx+1,1:Ny+1] ,100)
    pl.title(' $B_y(x, y)$')
    pl.xlabel('$x$ ')
    pl.ylabel('$y$')
    pl.colorbar()
    #pl.ylim([0,1])
    #pl.xlim([0,1])        
    #pl.grid(True)
    #pl.axes().set_aspect('equal', 'datalim')
    pl.savefig('magBy/point_mass' + '%04d'%time_index + '.png')
    pl.clf()    
    pl.close()
    
    pl.figure(figsize = (13,10) )
    pl.contourf(x_plot, y_plot ,div_B[1:Nx+1,1:Ny+1],100)
    pl.title( r'$\nabla \cdot \mathbf{B}$')
    pl.xlabel('$x$ ')
    pl.ylabel('$y$')
    pl.colorbar()
    #pl.ylim([0,1])
    #pl.xlim([0,1])        
    #pl.grid(True)
    #pl.axes().set_aspect('equal', 'datalim')
    pl.savefig('div/point_mass' + '%04d'%tim + '.png')
    #pl.show()
    pl.clf()   
    pl.close()
    
    pl.figure(figsize = (13,10) )
    pl.contourf(x_plot, y_plot ,Bx[1:Nx+1,1:Ny+1],100)
    pl.title(' $B_x(x, y)$')
    pl.xlabel('$x$ ')
    pl.ylabel('$y$')
    pl.colorbar()
    #pl.ylim([0,1])
    #pl.xlim([0,1])        
    #pl.grid(True)
    #pl.axes().set_aspect('equal', 'datalim')
    pl.savefig('magBx/point_mass' + '%04d'%time_index + '.png')
    pl.clf()         
    pl.close()
    

 
