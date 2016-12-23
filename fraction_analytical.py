from __future__ import division
import numpy as np
import pylab as pl
from scipy.integrate import nquad
from sympy import integrate, Symbol, exp


# Set plot parameters to make beautiful plots
pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20  
pl.rcParams['font.sans-serif'] = 'serif'
#pl.rcParams['text.usetex']     = True
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

l=1



    
def frac_above(t):
    v = Symbol('v')
    k=1
    m=1
    T=1
    f = np.sqrt(m/(2*np.pi*k*T))*exp(-m*v**2/(2*k*T))
    val = 2*integrate(f,(v,(l/(2*t)),np.inf))
    return val

def temp(t):
    if(t==0):
        return 0
    val = 2*frac_above(t)+(1-frac_above(t))*1.5
    return val


pl.figure()
pl.ylim(1.4,2.2)
for t in range(90):
    
    pl.plot(t,temp(0.1*t), 'o',color='blue', markersize=5, alpha = 0.4)
    pl.title('Plot')
    print(t)


pl.show()
pl.clf()
    
