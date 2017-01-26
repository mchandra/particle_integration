import numpy as np
import pylab as pl
import numpy.linalg as la
from interpolator import bilinear_interpolate
from fields import *

# Plotting parameters
pl.rcParams['figure.figsize']     = 12, 7.5
pl.rcParams['lines.linewidth']    = 1.5
pl.rcParams['font.family']        = 'serif'
pl.rcParams['font.weight']        = 'bold'
pl.rcParams['font.size']          = 20
pl.rcParams['font.sans-serif']    = 'serif'
pl.rcParams['text.usetex']        = True
pl.rcParams['axes.linewidth']     = 1.5
pl.rcParams['axes.titlesize']     = 'medium'
pl.rcParams['axes.labelsize']     = 'medium'
pl.rcParams['text.latex.unicode'] = True
pl.rcParams['xtick.major.size']   = 8
pl.rcParams['xtick.minor.size']   = 4
pl.rcParams['xtick.major.pad']    = 8
pl.rcParams['xtick.minor.pad']    = 8
pl.rcParams['xtick.color']        = 'k'
pl.rcParams['xtick.labelsize']    = 'medium'
pl.rcParams['xtick.direction']    = 'in'
pl.rcParams['image.aspect']       = 'equal'
pl.rcParams['ytick.major.size']   = 8
pl.rcParams['ytick.minor.size']   = 4
pl.rcParams['ytick.major.pad']    = 8
pl.rcParams['ytick.minor.pad']    = 8
pl.rcParams['ytick.color']        = 'k'
pl.rcParams['ytick.labelsize']    = 'medium'
pl.rcParams['ytick.direction']    = 'in'

"""user defined parameters"""

spread = 0.1

""" Number of ghost cells at each end (Keep it at 1)"""

ghost_cells = 1

""" Speed of light"""

c = 1

""" Time for the simulation"""

time_in_seconds = 2

""" Size of domain in x and y directions"""
Lx = 1
Ly = 1

"""Currents for the current fdtd code all irrelevant"""
Jx = 0
Jy = 0
Jz = 0

""" User defined functions for convenience """


def sumsum(a):

  return sum(sum(abs(a)))

def gauss2D(x,y):

  return np.exp(-( (x - Lx/2)**2 +(y - Ly/2)**2  )/(2*spread**2))

def gauss1D(x):

  return np.exp(-( (x - 0.5)**2 )/(2*spread**2))

def initial_fields(x, y):

  function_value = np.sin(2 * np.pi * x * y) * np.cos(2 * np.pi * x * y)

  return function_value
