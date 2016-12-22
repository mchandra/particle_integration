import numpy as np

""" Setting number of particles and other parameters"""

no_of_particles = 10000
x_divisions=32
y_divisions=1
z_divisions=1

""" Setting velocities according to maxwellian distribution """
k=1.0
m=1.0
T=1.5
T_walls=2
const=np.sqrt((k*T)/(m))
const2=np.sqrt((k*T_walls)/(m))
 
left_boundary = 0
right_boundary = 1
length_of_box_x         = right_boundary - left_boundary

bottom_boundary = 0
top_boundary = 1
length_of_box_y           = top_boundary - bottom_boundary

back_boundary = 0
front_boundary = 1
length_of_box_z           = front_boundary - back_boundary
