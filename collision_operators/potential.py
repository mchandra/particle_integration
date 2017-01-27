from simulation_parameters import *
from modules import *

def potential(x):
  potential= (potential_amplitude/2) * ( -1 * np.tanh(potential_steepness*x) + 1)
  return(potential)

# This function returns the value of potential gradient
# Here 'x' denotes the location where gradient is to be computed
# Here 'a' is the parameter which controls the steepness of the potential function

def potential_gradient(x, order):

  # We shall use the finite difference method to obtain the gradient of potential
  # There are several schemes of finite differencing with different orders of accuracy
  # Described below are the schemes used
  # In the below description O(deltaX) is the big-O notation 
  # In all the schemes, we have taken deltaX = 1e-10

  # In an order one scheme. The gradient of a quantity f w.r.t may be written as:
  # df/dx = ((f_{i+1} - f{i})/deltaX) + O(deltaX) 
  order1 = 1e10*(  potential(x+1e-10) - potential(x)     )

  # In an order one scheme. The gradient of a quantity f w.r.t may be written as:
  # df/dx = ((f_{i+1} - f{i-1})/2*deltaX) + O(deltaX^2)
  order2 = 5e9 *(  potential(x+1e-10)  - potential(x-1e-10))
  
  # In an order one scheme. The gradient of a quantity f w.r.t may be written as:
  # df/dx = (8*(f_{i+1} - f{i-1} + f_{i-2} - f_{i+2})/12*deltaX) + O(deltaX^4)
  order4 = (8  *(  potential(x+1e-10)  - potential(x-1e-10)) \
                 + potential(x-2e-10)  - potential(x+2e-10)  \
           ) /(12e-10)

  if   (order==4):
    return(order4)

  elif (order==2):
    return(order2)

  elif (order==1):
    return(order1)

  else:
    print("Error, order passed as argument not defined!")
    exit()

# This function calculates and returns the potential energy of the entire system
# This is done by summing over all the potentials of the particles in the system

def calc_potential_energy(sol):

  x = sol[0:no_of_particles].copy()               
  y = sol[no_of_particles:2*no_of_particles].copy() 

  # We shall use the fact that a * np.ones(a.size,a.size),generates a matrix
  # Of size a.size,a.size where each of the row of the matrix is the vector a
  # Thus it'll be of the form : array([[a],[a],[a],[a]......[a]])

  # Thus the following steps will generate a matrix, in which i-th row contains 
  # The difference in the x-coordinates of all particles - x-coordinate of i-th particle
  x = x * np.ones((no_of_particles,no_of_particles),dtype=np.float)
  x = x - np.transpose(x)

  # The similar transformation is performed for all y-coordinates
  y = y * np.ones((no_of_particles,no_of_particles),dtype=np.float)
  y = y - np.transpose(y)

  # dist is [N x N]. dis[i, j] is the distance between particle i and particle j.
  dist = np.sqrt(x**2+y**2)

  # potential(a,dist) will return a [N X N] matrix, where
  # potential[i,j] is the pair potential between particles i and j
  # Summing over all potentials in the system will give us the total potential energy
  # Of the system at a particular time-step
  potential_energy = 0.5*(np.sum(potential(dist)-(potential_amplitude/2)*np.identity(no_of_particles)))
  return(potential_energy)
