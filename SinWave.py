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

# dExdt = -dBydz
# dBydt = -dExdz

spread = 0.1


def gauss(zvar,time):
    
    return np.exp(  -  (zvar-0.5 )**2/(2*spread**2))  #(0.5+time*dt) - int(0.5+time*dt)

ghostcells = 1

def error(Nz):
    absSumErrorEx = 0
    absSumErrorBy = 0
    z = np.linspace(0,1,Nz+1)
    z_plot = z[0:Nz]
    Ex = np.zeros((Nz+2*ghostcells),dtype = np.float)
    By = np.zeros((Nz+2*ghostcells), dtype = np.float)

    div_B = np.zeros(Nz, dtype = np.float)
    c = 1


    dz = np.float(1/(Nz))
    dt = np.float(dz/(4*c))
    

    max_iterations = np.int(2/(dt))

    dt_by_dz = dt/(dz)
    for time_index in range(max_iterations):
        print('time = ', time_index,'grid size = ',Nz)
        if(time_index==0):
            for i in range(Nz):
                Ex[i+ghostcells] = np.sin(-2*np.pi*z[i])
                By[i+ghostcells] = np.sin(-2*np.pi*z[i])

        for i in range(ghostcells,Nz+ghostcells):

            Ex[i] = Ex[i] - (dt_by_dz*(By[i]-By[i-1]))


        
        Ex[0] = Ex[Nz]
        Ex[Nz+ghostcells] = Ex[ghostcells]
        
        
        
        
        
        for i in range(ghostcells,Nz+ghostcells):

            By[i] = By[i] - (dt_by_dz*(Ex[i+1] - Ex[i]))
        
        By[0] = By[Nz]
        By[Nz+ghostcells] = By[ghostcells]

        pl.plot(z_plot,Ex[1:Nz+1],'--',color = 'black',lw =3,label = 'Numerical ' )
        pl.legend()
        pl.plot(z_plot,np.sin(2*np.pi*(time_index*dt-z_plot)),color = 'green',lw =1,label = 'Analytical ' )
        pl.legend()
        pl.title(' $N  = ' + str(Nz)+'$')
        pl.xlabel('$z$ ')
        pl.ylabel('$E_x$')
        pl.ylim(-1.2,1.2)
        pl.legend()            
        pl.savefig('images/point_mass' + '%04d'%time_index + '.png')
        pl.clf()


        pl.plot(z_plot,Ex[1:Nz+1],'--',color = 'black',lw =3,label = 'Numerical ' )
        pl.legend()
        pl.plot(z_plot,np.sin(2*np.pi*(time_index*dt-z_plot)),color = 'green',lw =1,label = 'Analytical ' )
        pl.legend()
        pl.title(' $N  = ' + str(Nz)+'$')
        pl.xlabel('$z$ ')
        pl.ylabel('$B_y$')
        pl.ylim(-1.2,1.2)
        pl.legend()            
        pl.savefig('mag/point_mass' + '%04d'%time_index + '.png')
        pl.clf()

        if(  time_index==max_iterations-3  ):

            pl.plot(z_plot,Ex[1:Nz+1],'--',color = 'black',lw =3,label = 'Numerical ' )
            pl.legend()
            pl.plot(z_plot,np.sin(2*np.pi*(time_index*dt-z_plot)),color = 'green',lw =1,label = 'Analytical ' )
            pl.ylim(-1.2,1.2)
            pl.legend()            
            pl.show()
            pl.clf()

            absSumErrorEx = sum( abs(  Ex[ghostcells:Nz+ghostcells] - np.sin(2*np.pi*( (time_index+1)*dt-z_plot)  )  )   ) / (Nz)
            absSumErrorBy = sum( abs(  By[ghostcells:Nz+ghostcells] - np.sin(2*np.pi*( (time_index+1)*dt-z_plot) )  )   ) / (Nz)
        
            #print('Grid Points Taken =', Nz, ' Error in Electric Field= ', absSumErrorEx*Nz,' Error in Magnetic Field', absSumErrorBy*Nz)         
            return absSumErrorEx,absSumErrorBy


N = np.array( [ 128 ] )
ErrorNEx = np.zeros(len(N),dtype = np.float)
ErrorNBy = np.zeros(len(N),dtype = np.float)
for i in range(len(N)):
    ErrorNEx[i],ErrorNBy[i] = error(N[i])


pl.loglog(N,ErrorNEx,'-o',lw =3,label = '$E_x$ ' )
pl.legend()
pl.loglog(N,ErrorNBy,'-o',lw =3,label = '$B_y$ ' )
pl.legend()
pl.loglog(N,15*(N**-1.999),'--',color = 'black',lw = 2,label = ' $O(N^{-2})$. ')
pl.legend()
pl.title('$\mathrm{Convergence\; plot}$ ')
pl.xlabel('$\mathrm{N}$.')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
pl.show()
pl.clf()
    
