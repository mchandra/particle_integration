import pylab as pl
import numpy as np
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

""" Extremely basic plot file. This is to provide the reader an idea as to how 
data can be plotted. The needs specific to the user may be made here, like the 
time step for the simulation data which needs to be plotted. All files which 
were written in run.py can be accessed here, and can be plotted by changing the
options below"""

h5f = h5py.File('../data_files/timestepped_data/solution_500.h5', 'r')
pressure = h5f['pressure'][:]
time = h5f['time'][:]
h5f.close()

pl.plot(time, pressure)
pl.title('$\mathrm{Pressure}$ $\mathrm{Variation}$')
pl.xlabel('$\mathrm{Time}$')
pl.ylabel('$\mathrm{Pressure}$')
pl.savefig('plot.png')