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

a = np.array([-1.02522908456, -2.24048808396, -3.80990487832, -5.03251173974, -7.16147822735, -10.1289073704, -11.6426443717, -14.5094070456])
b = np.arange(1,9)
c = np.array([-22.5235889241,-59.5663055852,-110.633376964,-175.33967368,-253.360189853,-344.407640002,-448.220071002,-564.553187593])
c = c/32

pl.plot(b, a,label = 'Numerical')
pl.plot(b, c,label = 'Analytical')
pl.legend()
pl.title('$q_x$ $\mathrm{Versus}$ $\Delta T$')
pl.xlabel('$\Delta T$')
pl.ylabel('$\mathrm{HeatFlux}$')
pl.savefig('plot.png')