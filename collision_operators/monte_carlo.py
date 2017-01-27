from modules import *
from collision_parameters import *

def collision_operator(sol):
  x_zones     = np.zeros(sol.size,dtype=np.float)
  y_zones     = np.zeros(sol.size,dtype=np.float)

  count       = np.zeros((x_divisions,y_divisions),dtype=np.float)
  temperature = np.zeros((x_divisions,y_divisions),dtype=np.float)

  x_coordinates = sol[0:no_of_particles]
  y_coordinates = sol[no_of_particles:2*no_of_particles]                 
  velocity_x    = sol[2*no_of_particles:3*no_of_particles]
  velocity_y    = sol[3*no_of_particles:4*no_of_particles]

  for m in range(x_divisions):
    temporary1 = np.where(xCoords>x[m])
    temporary2 = np.where(xCoords<=x[m+1])
    indicesx = np.nonzero(np.in1d(temporary1, temporary2))[0]
    for n in range(y_divisions):
      temporary1 = np.where(yCoords[indicesx]>=y[n])
      temporary2 = np.where(yCoords[indicesx]<=y[n+1])
      indicesxy = np.nonzero(np.in1d(temporary1, temporary2))[0]
      collision[m][n] = indicesxy.tolist()
      count[m][n] = len(collision[m][n])
      temp[m][n] = temp[m][n] + 0.5*np.sum(velocityX[collision[m][n]]**2 + velocityY[collision[m][n]]**2)
    
  temp=temp/count
  tempglobal = temp
  for m in range(x_divisions):
    for n in range(y_divisions):
      x1 = np.random.rand(count[m][n])
      x2 = np.random.rand(count[m][n])
      sol[np.asarray(collision[m][n])+2*no_of_particles] = np.sqrt(2*temp[m][n])*erfinv(2*x1-1) 
      sol[np.asarray(collision[m][n])+3*no_of_particles] = np.sqrt(2*temp[m][n])*erfinv(2*x2-1) 

  return(sol)