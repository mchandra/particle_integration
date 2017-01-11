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


dx = np.float(1/(Nx))
dy = np.float(1/(Ny))
dt = np.float(dx/(2*c))


max_iterations = np.int(2/(dt))

dt_by_dx = dt/(dx)
dt_by_dy = dt/(dy)

for time_index in range(max_iterations):
    print('Time Index = ', time_index)
    if(time_index==0):
        for i in range(Nx):
            Ez[i+ghostcells][:] = np.exp(-(x[i]-0.5)**2/(2*spread**2))
            Bx[i+ghostcells][:] = 0
            By[i+ghostcells][:] = np.exp(-(x[i]-0.5)**2/(2*spread**2))

    for i in range(ghostcells, Nx + ghostcells):
        for j in range(ghostcells, Ny + ghostcells):
            
            Ez[i, j] = Ez[i, j] + (dt_by_dx*(By[i, j]-By[i-1, j])) -  (dt_by_dy*(Bx[i, j]-Bx[i,  j - 1 ]))


    
    Ez[0, :] = Ez[Nx, :].copy()
    Ez[: , 0] = Ez[:, Ny].copy()
    
    Ez[Nx+ghostcells, :] = Ez[ghostcells, :].copy()
    Ez[ :, Ny+ghostcells] = Ez[:, ghostcells].copy()
    
    
    
    
    for i in range(ghostcells,Nx+ghostcells):
        for j in range(ghostcells, Ny + ghostcells):
            
            Bx[i, j] = Bx[i, j] - (dt_by_dy*(Ez[i, j+1] - Ez[i, j]))
            
            By[i, j] = By[i, j] + (dt_by_dx*(Ez[i+1, j] - Ez[i, j]))
    
    
    
    Bx[0, :] = Bx[Nx, :].copy()
    Bx[:, 0] = Bx[:, Ny].copy()
    Bx[Nx+ghostcells, :] = Bx[ghostcells,:].copy()
    Bx[:, Ny+ghostcells] = Bx[:, ghostcells].copy()
    
    By[0, :] = By[Nx, :].copy()
    By[:, 0] = By[:, Ny].copy()
    By[Nx+ghostcells, :] = By[ghostcells,:].copy()
    By[:, Ny+ghostcells] = By[:, ghostcells].copy()
    
    for i in range(ghostcells,Nx+ghostcells):
        for j in range(ghostcells, Ny + ghostcells):
            
            div_B[i, j] = (Bx[i, j]-Bx[i-1,j])/(dx) +  (By[i, j]-By[i, j- 1])/(dy)

    div_B[0, :] = div_B[Nx, :].copy()
    div_B[:, 0] = div_B[:, Ny].copy()
    div_B[Nx+ghostcells, :] = div_B[ghostcells,:].copy()
    div_B[:, Ny+ghostcells] = div_B[:, ghostcells].copy()
    
    
    
    pl.contourf(Ez[1:Nx+1,1:Ny+1] )
    pl.title(' $E_x$ Movie with $N  = ' + str(Nx)+'$')
    pl.xlabel('$x$ ')
    pl.ylabel('$y$')
    pl.colorbar()
    pl.savefig('images/point_mass' + '%04d'%time_index + '.png')
    pl.clf()            
    
    
    pl.contourf(div_B[1:Nx+1,1:Ny+1])
    pl.title(' Divergence  Movie with $N  = ' + str(Nx)+'$')
    pl.xlabel('$x$ ')
    pl.ylabel
    pl.colorbar()
    pl.savefig('div/point_mass' + '%04d'%time_index + '.png')
    pl.clf()                
    
    
    #if(time_index==max_iterations-1):


        #absSumErrorEx = sum( abs(  Ex[ghostcells:Nz+ghostcells] - gauss(z_plot,time_index)  )   ) / (Nz)
        #absSumErrorBy = sum( abs(  By[ghostcells:Nz+ghostcells] - gauss(z_plot,time_index)  )   ) / (Nz)
#print('Grid Points Taken =', Nz, ' Error in Electric Field= ', absSumErrorEx,' Error in Magnetic Field', absSumErrorBy)         
#return absSumErrorEx,absSumErrorBy


#N = np.array( [ 32, 64, 128, 256, 512, 1024 ] )
#ErrorNEx = np.zeros(len(N),dtype = np.float)
#ErrorNBy = np.zeros(len(N),dtype = np.float)
#for i in range(len(N)):
    #ErrorNEx[i],ErrorNBy[i] = error(N[i])


#pl.loglog(N,ErrorNEx,'-o',lw =3,label = '$E_x$ ' )
#pl.legend()
#pl.loglog(N,ErrorNBy,'-o',lw =3,label = '$B_y$ ' )
#pl.legend()
#pl.loglog(N,15*(N**-1.999),'--',color = 'black',lw = 2,label = ' $O(N^{-2})$. ')
#pl.legend()
#pl.title('$\mathrm{Convergence\; plot}$ ')
#pl.xlabel('$\mathrm{N}$.')
#pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
#pl.show()
#pl.clf()
    
