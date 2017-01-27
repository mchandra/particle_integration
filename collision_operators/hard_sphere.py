from modules import *
from collision_parameters import *

def collision_operator(sol):
  if(simulation_dimension == 3):
    print("HS Scattering in 3D still in development! Please change to 2D mode")
  for i in range(no_of_particles):
    x_coordinates     = sol[(i+1):no_of_particles]        
    y_coordinates     = sol[(i+1+no_of_particles):2*no_of_particles]
    velocity_x        = sol[i+2*no_of_particles]*np.ones(x_coordinates.size)
    velocity_y        = sol[i+3*no_of_particles]*np.ones(y_coordinates.size)
    velocity_x_others = sol[(i+1+2*no_of_particles):3*no_of_particles]
    velocity_y_others = sol[(i+1+3*no_of_particles):4*no_of_particles]
    
    x_coordinates    = x_coordinates - sol[i]
    y_coordinates    = y_coordinates - sol[i+no_of_particles]
    dist             = np.sqrt(x_coordinates**2+y_coordinates**2) # Distance Vector, which will also be used in normalization
    x_coordinates    = x_coordinates/dist
    y_coordinates    = y_coordinates/dist

    test_collision = dist<0.01 
    test_collision = test_collision + 0
    indices        = np.nonzero(test_collision)
    
    if(np.sum(test_collision)!=0):
      p = (velocity_x*x_coordinates        + velocity_y*y_coordinates - \
           velocity_x_others*x_coordinates - velocity_y_others*y_coordinates\
          )*test_collision
      
      velocity_x        = velocity_x*test_collision - p*x_coordinates
      velocity_x_others = velocity_x_others + p*x_coordinates
      velocity_y        = velocity_y*test_collision - p*y_coordinates
      velocity_y_others = velocity_y_others + p*y_coordinates

      index = np.random.randint(0,(indices[0].size))
      
      sol[i+1+indices[0][index]+2*no_of_particles] = velocity_x_others[indices[0][index]]
      sol[i+1+indices[0][index]+3*no_of_particles] = velocity_y_others[indices[0][index]]
      sol[i+2*no_of_particles]                     = velocity_x[indices[0][index]]
      sol[i+3*no_of_particles]                     = velocity_y[indices[0][index]]

    return(sol)