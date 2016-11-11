import numpy as np
import pylab as pl

#%%

%matplotlib inline

""" Set plot parameters to make beautiful plots """
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['axes.linewidth']  = 1.5
pl.rcParams['axes.titlesize']  = 'medium'
pl.rcParams['axes.labelsize']  = 'medium'
pl.rcParams['text.usetex']     = True
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
#%%

no_of_particles = np.array([1000,10000,100000,1000000,10000000])
error_L1        = np.array([0.1445,0.0515,0.01426,0.005068,0.0009538])

#%%
pl.figure()
pl.xlabel('$\mathrm{N}$')
pl.ylabel('$L_1\;\mathrm{error}$')
pl.loglog(no_of_particles, error_L1, 'o-')
pl.loglog(no_of_particles, no_of_particles**-0.5, '--', color='black')
pl.xlim([100, 1e8])

