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

def error(Nz):
    absSumError = 0
    z = np.linspace(0,1,Nz+1)
    z_plot = z[0:Nz]
    Ex = np.zeros((Nz),dtype = np.float)
    By = np.zeros((Nz), dtype = np.float)

    div_B = np.zeros(Nz, dtype = np.float)
    c = 1


    dz = np.float(1/(Nz))
    dt = np.float(dz/(2*c))
    #print(dt)
    

    max_iterations = np.int(2/(dt))

    dt_by_dz = dt/(dz)
    for time_index in range(max_iterations):
        
        if(time_index==0):
            for i in range(Nz):
                Ex[i] = np.exp(-(z[i]-0.5)**2/(2*spread**2))
                By[i] = np.exp(-(z[i]-0.5)**2/(2*spread**2))

        for i in range(Nz):
            if(i==0):
                Ex[i] = Ex[i] - (dt_by_dz*(By[i]-By[i-1+Nz]))           
            else:
                Ex[i] = Ex[i] - (dt_by_dz*(By[i]-By[i-1]))


        for i in range(Nz):
            
            if(i==Nz-1):
                By[i] = By[i] - (dt_by_dz*(Ex[i+1-Nz] - Ex[i]))
            
            else:
                By[i] = By[i] - (dt_by_dz*(Ex[i+1] - Ex[i]))
    

        if(time_index==max_iterations-1):
            #pl.plot(z_plot,Ex,'--g',lw =3,label = 'Numerical ' )
            #pl.legend()
            #pl.plot(z_plot,gauss(z_plot,time_index),'b',lw = 2,label = ' Analytical ')
            #pl.legend()
            #pl.title('Numerical vs Analytical after 400 timesteps')
            #pl.xlabel('Z axis')
            #pl.ylabel('Ex : Electric Field in x direction varying along Z axis')
            #pl.ylim(-1,1)
            #pl.show()
            #pl.clf()            
            absSumError = sum( abs(  Ex - gauss(z_plot,time_index)  )   ) / (Nz)
    
    #print(absSumError)            
    return absSumError



N = np.array( [32,64,128,256,512, 1024 , 2048 , 4096 ] )
ErrorN = np.zeros(len(N),dtype = np.float)
for i in range(len(N)):
    ErrorN[i] = error(N[i])


pl.loglog(N,ErrorN,'--g',lw =3,label = '$\mathrm{L_1 Norm}$ ' )
pl.legend()
pl.loglog(N,1000*(N**-1.999),'b',lw = 2,label = ' Analytical ')
pl.legend()
pl.title('Error vs Grid Points after 400 timesteps')
pl.xlabel('Number of grid points')
pl.ylabel('$\mathrm{L_1 Norm}$')
pl.show()
pl.clf()
    
