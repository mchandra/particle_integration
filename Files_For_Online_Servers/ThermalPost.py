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
heatflux = h5f['heatflux'][:]
pressure = h5f['pressure'][:]
time = h5f['time'][:]
h5f.close()

pl.title('$\mathrm{Pressure/HeatFlux}$')
pl.xlabel('$t$')
pl.ylabel('$<v^2>$ $\mathrm{or}$ $<v^2v_x>$')
pl.plot(time[0:4998],heatflux[0:4998],'b',label='$\mathrm{HeatFlux}$')
pl.plot(time[0:4998],pressure[0:4998],'r',label='$\mathrm{Pressure}$')
pl.legend()
pl.savefig('plot.png')

