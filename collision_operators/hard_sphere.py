import numpy as np
from scipy.special import erfinv
import h5py
import params

"""
This model uses the HS model for scattering the particles.
By this model the particles are effectively scattered in the 
same way as to how one models collisions amongst billiard 
balls. By this model we shall condsider only 2-body collisions.
Multi-body collisions are treated as 2 body collisions, with the 
choice of the 2 bodies being random.
"""

"""Here we shall assign values as set in params"""

no_of_particles      = params.no_of_particles
simulation_dimension = params.simulation_dimension
restart_simulation   = params.restart_simulation
choice_integrator    = params.choice_integrator
collision_operator   = params.collision_operator

if(collision_operator == "hardsphere"):
  scattering_distance = params.scattering_distance

elif(collision_operator == "potential-based"):
  potential_steepness     = params.potential_steepness
  potential_amplitude     = params.potential_amplitude
  order_finite_difference = params.order_finite_difference

elif(collision_operator == "montecarlo"):
  x_zones            = params.x_zones
  y_zones            = params.y_zones
  scattered_fraction = params.scattered_fraction

mass_particle      = params.mass_particle
boltzmann_constant = params.boltzmann_constant
T_initial          = params.T_initial
wall_condition_x   = params.wall_condition_x
wall_condition_y   = params.wall_condition_y
wall_condition_z   = params.wall_condition_z

if(wall_condition_x == "thermal"):
  T_left_wall  = params.T_left_wall
  T_right_wall = params.T_right_wall

if(wall_condition_y == "thermal"):
  T_top_wall = params.T_top_wall
  T_bot_wall = params.T_bot_wall

if(wall_condition_z == "thermal"):
  T_front_wall = params.T_front_wall
  T_back_wall  = params.T_back_wall

left_boundary    = params.left_boundary
right_boundary   = params.right_boundary
length_box_x     = params.length_box_x

bottom_boundary  = params.bottom_boundary
top_boundary     = params.top_boundary
length_box_y     = params.length_box_y

back_boundary    = params.back_boundary
front_boundary   = params.front_boundary
length_box_z     = params.length_box_z

#Here we complete import of all the variable from the parameters file

def collision_operator(sol):

  if(simulation_dimension == 3):
    print("HS Scattering in 3D still in development! Please change to 2D mode")

  for i in range(no_of_particles):
    j=sol[(i+1):no_of_particles]        
    k=sol[(i+1+no_of_particles):2*no_of_particles]
    velx=sol[i+2*no_of_particles]*np.ones(j.size)
    vely=sol[i+3*no_of_particles]*np.ones(j.size)
    velx_others=sol[(i+1+2*no_of_particles):3*no_of_particles]
    vely_others=sol[(i+1+3*no_of_particles):4*no_of_particles]
    x_particle=sol[i]
    y_particle=sol[i+no_of_particles]        
    j=ne.evaluate("j-x_particle")
    k=ne.evaluate("k-y_particle")
    dist=ne.evaluate("sqrt(j**2+k**2)")
    j=j/dist
    k=k/dist
    test_collision=ne.evaluate("dist<0.01") 
    test_collision=ne.evaluate("where(test_collision,1,0)")
    indices=np.nonzero(test_collision)
    
    if(np.sum(test_collision)!=0):
      p=(velx*j+vely*k-velx_others*j-vely_others*k)*test_collision
      velx=velx*test_collision-p*j
      velx_others=velx_others+p*j
      vely=vely*test_collision-p*k
      vely_others=vely_others+p*k 
      index=np.random.randint(0,(indices[0].size))
      sol[i+1+indices[0][index]+2*no_of_particles] = velx_others[indices[0][index]]

      sol[i+1+indices[0][index]+3*no_of_particles] = vely_others[indices[0][index]]

      sol[i+2*no_of_particles]=velx[indices[0][index]]
      sol[i+3*no_of_particles]=vely[indices[0][index]]

  return(sol)