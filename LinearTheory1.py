import numpy as np
import pylab as pl
from scipy.integrate import odeint
import h5py

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

m = 1
K = 1
T = 1
# Fourier Mode
# k = 2 * np.pi
k = 1

def f_0(v):
    return np.sqrt(m/(2*np.pi*K*T))*np.exp(-m*v**2/(2*K*T))

def diff_f_0_v(v):
    return np.sqrt(m/(2*np.pi*K*T))*np.exp(-m*v**2/(2*K*T)) * ( -m * v / (K * T))

def diff_Df(Y,t, v, int_Df_r, int_Df_i):
    f_r, f_i = Y
    # print('The first dYdt Term is ', (k * v * f_i) - (m * int_Df_i * diff_f_0_v(v)/k ))
    # print('The second dYdt Term is ', -(k * v * f_r) + (m * int_Df_r * diff_f_0_v(v)/k ))
    dYdt =[\
             (k * v * f_i) - (m * int_Df_i * diff_f_0_v(v)/k ),\
            -(k * v * f_r) + (m * int_Df_r * diff_f_0_v(v)/k )\
          ]
    # print(dYdt)
    return dYdt


velocity_max =  +10
velocity_min =  -10

number_of_velocities_points = 100
velocity_x = np.linspace(velocity_min, velocity_max, number_of_velocities_points)
# velocity_x = np.array([1])


# D_f_initial = (0.01) * np.exp(1j * velocity_x)
D_f_initial = 0.01 * f_0(velocity_x)

final_time = 8
dt = 0.001
time = np.arange(0, final_time, dt)

initial_conditions = np.zeros((len(velocity_x), 2), dtype = np.float)

max_amplitude = np.zeros(len(time), dtype = np.float)

for time_index, t0 in enumerate(time):

    print("Computing for TimeIndex = ", time_index)
    t0 = time[time_index]
    if (time_index == time.size - 1):
        break
    t1 = time[time_index + 1]
    t = [t0, t1]

    int_Df_r = 0
    int_Df_i = 0
    diff_f_all = np.zeros(len(velocity_x), dtype=np.complex128)

    if(time_index == 0):

        initial_conditions[:, 0] = D_f_initial.real
        initial_conditions[:, 1] = D_f_initial.imag

        int_Df_r = np.sum(D_f_initial.real)
        int_Df_i = np.sum(D_f_initial.imag)

    else:

        initial_conditions[:, 0] = old.real
        initial_conditions[:, 1] = old.imag

        int_Df_r = np.sum(old.real)
        int_Df_i = np.sum(old.imag)

    # print('D_f_initial are',D_f_initial)
    # print('initial_conditions are',initial_conditions)
    # print('int_Df_r are', int_Df_r)
    # print('int_Df_i are', int_Df_i)

    # print('time is ', t)
    # print('initial_conditions is ',initial_conditions )
    for i in range(len(velocity_x)):

        sol = odeint(diff_Df, initial_conditions[i], t, args=(velocity_x[i],int_Df_r, int_Df_i))[1]
        # print('SOlution is ',sol)
        diff_f_all[i] += sol[0] + 1j * sol[1]
            # print('inside ')


    # h5f = h5py.File('data_files/solution' + str(time_index) + '.h5', 'w')
    # h5f.create_dataset('data_files/solution_dataset' + str(time_index), data=diff_f_all)
    # h5f.close()

    # print('Solution', D_f_initial)
    max_amplitude[time_index] = sum(diff_f_all.real)
    # print('Solution', max(diff_f_all.real))
    old = diff_f_all.copy()
    # zzz = input('Whats up')


pl.plot(time, max_amplitude)
pl.xlabel('$\mathrm{time}$')
pl.ylabel(r'$\rho\left(t\right)$')
pl.title('$\mathrm{Linear\;decay}$')
pl.show()
