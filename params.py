import numpy as np
import pylab as pl
import numpy.linalg as la
from interpolator import bilinear_interpolate
from fields import *
import h5py

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

""" Speed of light (in natural units)"""

c = 1

""" Time for the simulation (in natural units) """

time_in_seconds = 2

""" Size of domain in x and y directions (in natural units)"""

Lx = 1
Ly = 1

"""Currents for the current fdtd code all irrelevant"""

Jx = 0
Jy = 0
Jz = 0

"""Charge/Mass of the particle"""
# mass and charge of particle:

me = 1
charge = 1


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


# Charge Deposition for B0 splines (Have to vectorize)

def b0_one_d(x_i,x_grid,dx):
    charged_grid = np.zeros(len(x_grid), dtype = np.float)
    # len(x_i) = number of particles
    for particle_index in len(x_i):
        for grid_index in len(x_grid):
            if(abs((xi[particle_index] - x_grid[grid_index])/dx)<0.5):
                charged_grid[grid_index]+=1
    return charged_grid



def b0_two_d(xi, yi, x_grid,y_grid, dx):
    charged_grid = np.zeros((len(x_grid), len(y_grid)),dtype = np.float)
    for particle_index in len(xi):
        for x_index in len(x_grid):
            for y_index in len(y_grid):
                if(abs((xi[particle_index]-x_grid[x_index])/dx<0.5) and abs((yi[particle_index]-y_grid[y_index])/dy<0.5) ):
                    charged_grid[x_index, y_index] +=1

    return charged_grid

def b0_three_d(xi, yi, x_grid,y_grid, z_grid, dx, dy, dz):
    charged_grid = np.zeros((len(x_grid), len(y_grid)),dtype = np.float)
    for particle_index in len(xi):
        for x_index in len(x_grid):
            for y_index in len(y_grid):
                for z_index in len(z_grid):
                    if(abs((xi[particle_index]-x_grid[x_index])/dx<0.5) and abs((yi[particle_index]-y_grid[y_index])/dy<0.5) \
                        and abs((zi[particle_index]-z_grid[z_index])/dz<0.5) \
                      ):
                        charged_grid[x_index, y_index, z_index] +=1

    return charged_grid


# Charge Deposition for B1 splines (Under Progress)

def b0_two_d(xyi, x_grid, y_grid, z_grid, dx):
    charged_grid = np.zeros(len(x_grid), dtype = np.float)
    charged_indices = np.where(abs((x_grid-xyi[:,0])/dx)<0.5 and (y_grid-xyi[:,1])/dx)<0.5 and (z_grid-xyi[:,1])/dx)<0.5):
    charged_grid[charged_indices]+= charge


def b1(xi,x_grid,dx):
    if(abs((x-x0)/dx)<1):
        return np.float(1-abs((x-x0)/dx))
    else:
        return 0

#
# b0 = np.vectorize(b0, otypes = [np.float], excluded = ['x_grid'])
# b1 = np.vectorize(b1, otypes = [np.float], excluded = ['x_grid'])
