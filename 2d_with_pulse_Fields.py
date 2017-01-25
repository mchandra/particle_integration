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

x_divisions = 100
y_divisions = 100
max_iterations = 1000

c0 = 3.e8
dx = 0.01
dt = dx/(2.*c0)

cc = (c0*dt/dx)

Ez = np.zeros((x_divisions,y_divisions),dtype = np.float)

Hx = np.zeros((x_divisions,y_divisions),dtype = np.float)
Hy = np.zeros((x_divisions,y_divisions),dtype = np.float)

spread = 8
t0 = 20

for time_index in range(max_iterations):
    print('Time Index', time_index)
    for k2 in range(1,y_divisions-1):
        for k1 in range(1,x_divisions-1):
            
            Ez[k1][k2] =  Ez[k1][k2] +  cc*(Hy[k1][k2]-Hy[k1-1][k2])+cc*(Hx[k1][k2-1]-Hx[k1][k2])                                   # Updating Ez
# Updating the Electric field
    
    
    Ez[int(x_divisions/2),int(y_divisions/2)]= np.exp(-0.5*((time_index-t0)/spread)**2)
# the pulse    
    
    for k2 in range(1,y_divisions-2):
        for k1 in range(1,x_divisions-2):
            Hx[k1][k2] = Hx[k1][k2] + cc*(Ez[k1][k2]-Ez[k1][k2+1])
            Hy[k1][k2] = Hy[k1][k2] + cc*(Ez[k1+1][k2]-Ez[k1][k2])
        
   # Updating the Magnetic Fields 
    
    
    pl.contourf(Ez,100)
    pl.savefig('images/point_mass' + '%04d'%time_index + '.png')
    pl.title('$\mathrm{Ez vs Time with Pulse}$ $\mathrm{Solution}$')
    pl.clf



