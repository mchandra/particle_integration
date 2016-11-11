import numpy as np
import pylab as pl
from scipy.integrate import odeint
import scipy.stats as stats
import matplotlib.pyplot as plt
import math


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

no_of_particles=[1000,10000,100000,1000000,10000000]
mean=[0.1445,0.0515,0.01426,0.005068,0.0009538]
variance=[0.0313806451613,0.00439857548387,0.000365383019355,4.05846874839e-05,1.35727896774e-06]
std_dev=[0.17714583021141178,0.06632175724353936,0.019114994620842513,0.006370611233144817,0.0011650231618907738]

pl.figure()
pl.title('Mean of Error vs No of Particles')
pl.xlabel('No of particles')
pl.ylabel('Mean')
pl.loglog(no_of_particles,mean)
pl.show()

pl.figure()
pl.title('Variance of Error vs No of Particles')
pl.xlabel('No of particles')
pl.ylabel('Variance')
pl.loglog(no_of_particles,mean)
pl.show()

pl.figure()
pl.title('Standard Deviation of Error vs No of Particles')
pl.xlabel('No of particles')
pl.ylabel('Standard Deviation')
pl.loglog(no_of_particles,mean)
pl.show()


