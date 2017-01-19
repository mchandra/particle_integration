import numpy as np
import h5py

# All quantities of length are mentioned in metres
noOfParticles  = 100
massParticle   = 5  
leftBoundary   = 0.
rightBoundary  = 1.
bottomBoundary = 0.
topBoundary    = 1.
lengthBoxX     = rightBoundary - leftBoundary 
lengthBoxY     = topBoundary   - bottomBoundary   

# Randomizing the initial conditions of the particles to spread them evenly through the box
initialPositionX = lengthBoxX * np.random.rand(noOfParticles)
initialPositionY = lengthBoxY * np.random.rand(noOfParticles)

# Assigning velocities so that a flat velocity distribution is obtained
initialVelocityX = np.random.rand(noOfParticles)
initialVelocityY = np.random.rand(noOfParticles)

# We shall randomize direction of velocity
initialVelocityX = initialVelocityX * np.random.choice([-1,1],size=noOfParticles)
initialVelocityY = initialVelocityY * np.random.choice([-1,1],size=noOfParticles)

# Combinining the initial conditions into a single vector
initialConditions = np.concatenate([initialPositionX,initialPositionY, \
                                    initialVelocityX,initialVelocityY, \
                                   ], axis = 0 \
                                  )

# All quantities of time are mentioned in seconds
boxCrossingTimeScale = (lengthBoxX/np.max(initialVelocityX))
finalTime            = 1.0 * boxCrossingTimeScale
dt                   = 0.0005 * boxCrossingTimeScale
time                 = np.arange(0, finalTime, dt)

# Declaring the arrays which will be used to store data

momentumX       = np.zeros(time.size)
momentumY       = np.zeros(time.size)

potentialEnergy = np.zeros(time.size)
kineticEnergy   = np.zeros(time.size)

# We define the pair potential that exists between any 2 particles in the medium
# Here 'x' denotes the distance between the two particles
# Here 'a' is the parameter which controls the steepness of the potential function

def potential(a, x):
  potential= 200 * ( -1 * np.tanh(a*x) + 1)
  return(potential)

# This function returns the value of potential gradient
# Here 'x' denotes the location where gradient is to be computed
# Here 'a' is the parameter which controls the steepness of the potential function

def potentialGradient(a, x, order):

  # We shall use the finite difference method to obtain the gradient of potential
  # There are several schemes of finite differencing with different orders of accuracy
  # Described below are the schemes used
  # In the below description O(deltaX) is the big-O notation 
  # In all the schemes, we have taken deltaX = 1e-10

  # In an order one scheme. The gradient of a quantity f w.r.t may be written as:
  # df/dx = ((f_{i+1} - f{i})/deltaX) + O(deltaX) 
  order1 = 1e10*(  potential(a, x+1e-10) - potential(a, x)     )

  # In an order one scheme. The gradient of a quantity f w.r.t may be written as:
  # df/dx = ((f_{i+1} - f{i-1})/2*deltaX) + O(deltaX^2)
  order2 = 5e9 *(  potential(a,x+1e-10)  - potential(a,x-1e-10))
  
  # In an order one scheme. The gradient of a quantity f w.r.t may be written as:
  # df/dx = (8*(f_{i+1} - f{i-1} + f_{i-2} - f_{i+2})/12*deltaX) + O(deltaX^4)
  order4 = (8  *(  potential(a,x+1e-10)  - potential(a,x-1e-10)) \
                 + potential(a,x-2e-10)  - potential(a,x+2e-10)  \
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

def calcPotentialEnergy(soln):

  x = soln[0:noOfParticles].copy()               
  y = soln[noOfParticles:2*noOfParticles].copy() 

  # We shall use the fact that a * np.ones(a.size,a.size),generates a matrix
  # Of size a.size,a.size where each of the row of the matrix is the vector a
  # Thus it'll be of the form : array([[a],[a],[a],[a]......[a]])

  # Thus the following steps will generate a matrix, in which i-th row contains 
  # The difference in the x-coordinates of all particles - x-coordinate of i-th particle
  x = x * np.ones((noOfParticles,noOfParticles),dtype=np.float)
  x = x - np.transpose(x)

  # The similar transformation is performed for all y-coordinates
  y = y * np.ones((noOfParticles,noOfParticles),dtype=np.float)
  y = y - np.transpose(y)

  # dist is [N x N]. dis[i, j] is the distance between particle i and particle j.
  dist = np.sqrt(x**2+y**2)

  # potential(a,dist) will return a [N X N] matrix, where
  # potential[i,j] is the pair potential between particles i and j
  # Summing over all potentials in the system will give us the total potential energy
  # Of the system at a particular time-step
  potentialEnergy = np.sum(potential(300,dist))
  return(potentialEnergy)


def Verlet(a,initialConditions,dt):

  x         = initialConditions[0:noOfParticles].copy()
  y         = initialConditions[noOfParticles:2*noOfParticles].copy()
  velocityX = initialConditions[2*noOfParticles:3*noOfParticles].copy()
  velocityY = initialConditions[3*noOfParticles:4*noOfParticles].copy()

  # Using the transformation which is used in calcPotentialEnergy()
  x = x * np.ones((noOfParticles,noOfParticles),dtype=np.float)
  x = x - np.transpose(x)

  y = y * np.ones((noOfParticles,noOfParticles),dtype=np.float)
  y = y - np.transpose(y)


  dist = np.sqrt(x**2+y**2)

  # Here vector is a 3D array, which is used to assign components of force
  # This is done by normalizing this vector to a unit vector, where
  # The array elements are used to describe the components of the unit vector
  # Of the line that joins any two particles
  # nvector[0] describes a 2D array, which has nvector[0][i,j],describe the 
  # x-component of the unit vector joining particle i and particle j
  # Similarly, nvector[1] describes a 2D array, which has nvector[1][i,j],describe the 
  # y-component of the unit vector joining particle i and particle j
  # np.nan_to_num is used to avoid 0/0 errors, where 0/0 gives 0 
  vector  = np.array([x,y])
  nvector = np.nan_to_num(vector/dist)

  # force is the force matrix, in which
  # force[i,j] denotes the force on particle i due to particle j
  force = potentialGradient(a,dist,4)

  # force X, and force Y are the 1D arrays in which
  # The i-th element describes the net-force on the i-th particle
  # In the x-direction and y-direction respectively
  # These are then used to correct the velocity components

  forceX    = np.sum(force * nvector[0],axis=1)
  velocityX = velocityX + 0.5*(forceX/massParticle)*dt

  forceY    = np.sum(force * nvector[1],axis=1)
  velocityY = velocityY + 0.5*(forceY/massParticle)*dt

  x = initialConditions[0:noOfParticles].copy()
  y = initialConditions[noOfParticles:2*noOfParticles].copy()

  xNew = x + velocityX*dt
  yNew = y + velocityY*dt

  x = xNew.copy()
  y = yNew.copy()
  
  #Employing the same transformation used earlier:

  x = x * np.ones((noOfParticles,noOfParticles),dtype=np.float)
  x = x - np.transpose(x)

  y = y * np.ones((noOfParticles,noOfParticles),dtype=np.float)
  y = y - np.transpose(y)

  dist    = np.sqrt(x**2+y**2)
  vector  = np.array([x,y])
  nvector = np.nan_to_num(vector/dist)

  forceX    = np.sum(force * nvector[0],axis=1)
  velocityX = velocityX + 0.5*(forceX/massParticle)*dt

  forceY    = np.sum(force * nvector[1],axis=1)
  velocityY = velocityY + 0.5*(forceY/massParticle)*dt

  nextStep=np.concatenate([xNew, yNew, velocityX, velocityY],axis=0)
  return(nextStep)

# sol is the vector that describes the current state of the system
# sol consists of the x-coordinate, y-coordinate and
# The x and y components of velocity for all the particles.
# old is used to store the sol vector of the previous time step
sol = np.zeros(4*noOfParticles,dtype=np.float)
old = np.zeros(4*noOfParticles,dtype=np.float)

# Now we shall proceed to evolve the system with time
for timeIndex,t0 in enumerate(time):
  
  print("Computing For Time Index = ",timeIndex)
  
  if(timeIndex == time.size-1):
    break
  
  if(timeIndex==0):
    initialConditions = initialConditions
  else:
    initialConditions = old
  
  sol = Verlet(300,initialConditions,dt)

  # The following section implements periodic B.C's
  # It works by recognizing where particles have crossed the walls
  # And then corrects those indices where the condition is observed to be true
  xCoords = sol[0:noOfParticles]
  yCoords = sol[noOfParticles:2*noOfParticles]

  wallYBot   = np.where(yCoords<bottomBoundary)
  wallYTop   = np.where(yCoords>topBoundary)
  wallXLeft  = np.where(xCoords<leftBoundary)
  wallXRight = np.where(xCoords>rightBoundary)

  yCoords[wallYBot[0]]   = yCoords[wallYBot[0]]   + lengthBoxY
  yCoords[wallYTop[0]]   = yCoords[wallYTop[0]]   - lengthBoxY
  xCoords[wallXLeft[0]]  = xCoords[wallXLeft[0]]  + lengthBoxX
  xCoords[wallXRight[0]] = xCoords[wallXRight[0]] - lengthBoxX

  sol[noOfParticles+wallYBot[0]] = yCoords[wallYBot[0]]
  sol[noOfParticles+wallYTop[0]] = yCoords[wallYTop[0]]
  sol[wallXLeft[0]]              = xCoords[wallXLeft[0]]
  sol[wallXRight[0]]             = xCoords[wallXRight[0]]

  old = sol

  # Computing total momentum and energies of the system
  totalMomentumX       = massParticle * np.sum(sol[(2*noOfParticles):(3*noOfParticles)])
  totalMomentumY       = massParticle * np.sum(sol[(3*noOfParticles):(4*noOfParticles)])
  totalKineticEnergy   = 0.5*massParticle*np.sum(sol[2*noOfParticles:3*noOfParticles]**2 +\
                                                 sol[3*noOfParticles:4*noOfParticles]**2
                                                )
  totalPotentialEnergy = calcPotentialEnergy(sol)

  # Assigning the values to an array, which will be used in post-processing
  momentumX[timeIndex]       = totalMomentumX
  momentumY[timeIndex]       = totalMomentumY
  kineticEnergy[timeIndex]   = totalKineticEnergy
  potentialEnergy[timeIndex] = totalPotentialEnergy
  
  # Writing the data to file every 1000 time steps
  # This data will then be post-processed to generate results
  if((timeIndex%1000)==0):
    h5f = h5py.File('solution_data/solution_'+str(timeIndex)+'.h5', 'w')
    h5f.create_dataset('sol',             data=sol)
    h5f.create_dataset('momentumX',       data=momentumX)
    h5f.create_dataset('momentumY',       data=momentumY)
    h5f.create_dataset('kineticEnergy',   data=kineticEnergy)
    h5f.create_dataset('potentialEnergy', data=potentialEnergy)
    h5f.create_dataset('time',            data=time)
    h5f.close()
