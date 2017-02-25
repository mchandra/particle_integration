# REFER TO THE LATEX DOCUMENT FOR DETAILS ON THE DIFFERENTIAL EQUATIONS TO BE INTEGRATED
import numpy as np
import pylab as pl
from scipy.integrate import odeint
import h5py


# In[ ]:

pl.rcParams['figure.figsize']  = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family']     = 'serif'
pl.rcParams['font.weight']     = 'bold'
pl.rcParams['font.size']       = 20
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex']     = True
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


# In[ ]:

# Setting the variables in the maxwell distribution
m = 20
K = 1
T = 1
# In[ ]:

# k for the mode in fourier space
k = 1


# In[ ]:

# The maxwell Boltzman function
def f_0(v):
    return np.sqrt(m/(2*np.pi*K*T))*np.exp(-m*v**2/(2*K*T))

# This the function which returns the derivative of the maxwell boltzmann equation
def diff_f_0_v(v):
    return np.sqrt(m/(2*np.pi*K*T))*np.exp(-m*v**2/(2*K*T)) * ( -m * v / (K * T))


# In[ ]:

# Assign the maxim and minimum velocity for the velocity grid
velocity_max =  +10
velocity_min =  -10

# Set the divisions for the velocity grid
number_of_velocities_points = 101
velocity_x = np.linspace(velocity_min, velocity_max, number_of_velocities_points)
dv = velocity_x[1] - velocity_x[0]


# In[ ]:

# Function that returns df_i/dt and df_r/dt used for odeint function
# See the latex document for more details on the differential equations
# This has been done to split the imaginary and real part of the ODE
def diff_delta_f(Y,t):
    f_r = Y[0:len(velocity_x)]  # Initial conditions for odeint
    f_i = Y[len(velocity_x): 2 * len(velocity_x)] 
    
    int_Df_i = np.sum(f_i) * (velocity_x[1]-velocity_x[0])
    int_Df_r = np.sum(f_r) * (velocity_x[1]-velocity_x[0])

    # This the derivate for f_r and f_i given in the latex document
    dYdt =np.concatenate([(k * velocity_x * f_i) - (int_Df_i * diff_f_0_v(velocity_x)/k ), \
                           -(k * velocity_x * f_r) + (int_Df_r * diff_f_0_v(velocity_x)/k )\
                         ], axis = 0)
    # This returns the derivative for the coupled set of ODE

    return dYdt

def diff_delta_f_Ex(Y,t):

    f_r = Y[0:len(velocity_x)]  # Initial conditions for odeint
    f_i = Y[len(velocity_x): 2 * len(velocity_x)]
    E_x_r = Y[2 * len(velocity_x)] 
    E_x_i = Y[2 * len(velocity_x) + 1]
    #print('Er', E_x_r)
    #print('Ei', E_x_i)
    #input('whats up')
    int_v_delta_f_dv_i = np.sum(f_i * velocity_x) * (dv)
    int_v_delta_f_dv_r = np.sum(f_r * velocity_x) * (dv)
    int_v_delta_f_dv = np.array([int_v_delta_f_dv_r, int_v_delta_f_dv_i ] )

    # This the derivate for f_r and f_i given in the latex document
    dYdt =np.concatenate([(    k * velocity_x * f_i) - (E_x_r * diff_f_0_v(velocity_x) ), \
                            - (k * velocity_x * f_r) - (E_x_i * diff_f_0_v(velocity_x) ), \
                                -1 * int_v_delta_f_dv\
                         ], axis = 0\
                        )
    # This returns the derivative for the coupled set of ODE
    
    return dYdt
    
# In[ ]:

# Set the initial conditions for delta f(v,t) here
delta_f_initial = np.zeros((2 * len(velocity_x)), dtype = np.float)
delta_f_initial[0: len(velocity_x)] = 0.01 * f_0(velocity_x)

delta_f_Ex_initial = np.zeros((2 * len(velocity_x)+2), dtype = np.float)
delta_f_Ex_initial[0 : len(velocity_x)] = 0.01 * f_0(velocity_x)
delta_f_Ex_initial[2 * len(velocity_x) + 1] = -1 * (1/k) * np.sum(delta_f_Ex_initial[0: len(velocity_x)] ) * dv

# In[ ]:

# Setting the parameters for time here
final_time = 10
dt = 0.001
time = np.arange(0, final_time, dt)


# In[ ]:

# Variable for temperorily storing the real and imaginary parts of delta f used for odeint
initial_conditions_delta_f = np.zeros((2 * len(velocity_x)), dtype = np.float)
old_delta_f = np.zeros((2 * len(velocity_x)), dtype = np.float)


initial_conditions_delta_f_Ex = np.zeros((2 * len(velocity_x) + 2), dtype = np.float)
old_delta_f_Ex = np.zeros((2 * len(velocity_x) + 2 ), dtype = np.float)
# Variable for storing delta rho

delta_rho1 = np.zeros(len(time), dtype = np.float)
delta_rho2 = np.zeros(len(time), dtype = np.float)
delta_f_temp = np.zeros(2 * len(velocity_x), dtype=np.float)
temperory_delta_f_Ex = np.zeros(2 * len(velocity_x) + 2, dtype=np.float)
# In[ ]:

for time_index, t0 in enumerate(time):
    if(time_index%1000==0):
        print("Computing for TimeIndex = ", time_index)
    t0 = time[time_index]
    if (time_index == time.size - 1):
        break
    t1 = time[time_index + 1]
    t = [t0, t1]

    # delta f is defined on the velocity grid


    # Initial conditions for the odeint 
    if(time_index == 0):
        # Initial conditions for the odeint for the 2 ODE's respectively for the first time step
        # First column for storing the real values of delta f and 2nd column for the imaginary values
        initial_conditions_delta_f                 = delta_f_initial.copy()
        initial_conditions_delta_f_Ex                 = delta_f_Ex_initial.copy()
        # Storing the integral sum of delta f dv used in odeint

    else:
        # Initial conditions for the odeint for the 2 ODE's respectively for all other time steps
        # First column for storing the real values of delta f and 2nd column for the imaginary values
        initial_conditions_delta_f= old_delta_f.copy()
        initial_conditions_delta_f_Ex= old_delta_f_Ex.copy()
        # Storing the integral sum of delta f dv used in odeint

    # Integrating delta f
    
    temperory_delta_f = odeint(diff_delta_f, initial_conditions_delta_f, t)[1]
    temperory_delta_f_Ex = odeint(diff_delta_f_Ex, initial_conditions_delta_f_Ex, t)[1]

    # Saving delta rho for current time_index
    delta_rho1[time_index] = ((sum(dv * temperory_delta_f[0: len(velocity_x)])))
    delta_rho2[time_index] = ((sum(dv * temperory_delta_f_Ex[0: len(velocity_x)])))
    
    # Saving the solution for to use it for the next time step
    old_delta_f = temperory_delta_f.copy()
    old_delta_f_Ex = temperory_delta_f_Ex.copy()
    #if(time_index%100==0):
        #print(temperory_delta_f[0: len(velocity_x)])
        #input('test')
# In[ ]:

# Plotting the required quantities here
pl.plot(time, np.log(abs(delta_rho1)), '-o' ,label = '$\mathrm{First\;Approach}$')
pl.plot(time, np.log(abs(delta_rho2)), '--' ,label = '$\mathrm{Second\;Approach}$')
pl.xlabel('$\mathrm{time}$')
pl.ylabel(r'$\delta \hat{\rho}\left(t\right)$')
pl.title('$\mathrm{Linear\;Landau\;damping}$')
pl.legend()
pl.show()


# In[ ]:



