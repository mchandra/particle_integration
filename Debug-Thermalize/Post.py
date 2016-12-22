import numpy as np
import pylab as pl
import h5py

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

h5f = h5py.File('post.h5', 'r')
heatflux_x = h5f['heatflux_x'][:]
heatflux_y = h5f['heatflux_y'][:]
heatflux_z = h5f['heatflux_z'][:]
pressure = h5f['pressure'][:]
time = h5f['time'][:]
h5f.close()

pl.title('$\mathrm{Pressure/HeatFlux}$')
pl.xlabel('$t$')
pl.plot(time[:-1],heatflux_x[:-1],'b',label='$<v^2v_x>$')
pl.plot(time[:-1],heatflux_y[:-1],'g',label='$<v^2v_y>$')
pl.plot(time[:-1],heatflux_z[:-1],'y',label='$<v^2v_z>$')
pl.plot(time[:-1],pressure[:-1],'r',label='$<v^2>$')
pl.legend(loc='center right')
pl.savefig('plot.png')

