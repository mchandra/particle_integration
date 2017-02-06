import params

simulation_dimension = params.simulation_dimension


if(simulation_dimension == 3):

  def collision_operator(xcoords, ycoords, zcoords, vel_x, vel_y, vel_z, dt):
   return(xcoords, ycoords, zcoords, vel_x, vel_y, vel_z)

if(simulation_dimension == 2):

  def collision_operator(xcoords, ycoords, vel_x, vel_y, dt):
   return(xcoords, ycoords, vel_x, vel_y)