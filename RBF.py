import numpy as np
import pylab as pl
import numpy.linalg as la

pl.rcParams['figure.figsize'] = 12, 7.5
pl.rcParams['lines.linewidth'] = 1.5
pl.rcParams['font.family'] = 'serif'
pl.rcParams['font.weight'] = 'bold'
pl.rcParams['font.size'] = 20
pl.rcParams['font.sans-serif'] = 'serif'
pl.rcParams['text.usetex'] = True
pl.rcParams['axes.linewidth'] = 1.5
pl.rcParams['axes.titlesize'] = 'medium'
pl.rcParams['axes.labelsize'] = 'medium'

pl.rcParams['xtick.major.size'] = 8
pl.rcParams['xtick.minor.size'] = 4
pl.rcParams['xtick.major.pad'] = 8
pl.rcParams['xtick.minor.pad'] = 8
pl.rcParams['xtick.color'] = 'k'
pl.rcParams['xtick.labelsize'] = 'medium'
pl.rcParams['xtick.direction'] = 'in'

pl.rcParams['ytick.major.size'] = 8
pl.rcParams['ytick.minor.size'] = 4
pl.rcParams['ytick.major.pad'] = 8
pl.rcParams['ytick.minor.pad'] = 8
pl.rcParams['ytick.color'] = 'k'
pl.rcParams['ytick.labelsize'] = 'medium'
pl.rcParams['ytick.direction'] = 'in'


def get_xi(nodal_point, local_domain_nodes):
    return np.float(nodal_point / (local_domain_nodes - 1))


def shape(localXi, local_array):
    # print('The Xi for the term is ', localXi)
    # print('Number of points provided =',len(local_array))
    local_domain_nodes = len(local_array)
    spread = np.float((1 / local_domain_nodes)+0.1)
    localXi_matrix = [[0] for i in range(local_domain_nodes)]
    u = np.matrix(local_array)

    for i in range(local_domain_nodes):
        localXi_matrix[i] = np.exp(-  (localXi - get_xi(i, local_domain_nodes)) ** 2 / (2 * spread ** 2))
    # print('local Xi vector matrix is ', localXi_matrix)
    c = np.matrix([[0] for i in range(local_domain_nodes)])
    A = np.zeros((local_domain_nodes, local_domain_nodes), dtype=np.float)

    for i in range(local_domain_nodes):
        for j in range(local_domain_nodes):
            A[i, j] = np.exp(- (get_xi(i, local_domain_nodes) - get_xi(j, local_domain_nodes)) ** 2 / (2 * spread ** 2))
    # print('A is  ', A)
    c = np.matrix(la.inv(A)) * u.transpose()
    # print('c is  is ', c)
    # print('The answer returned is  ', np.matrix(localXi_matrix) * np.matrix(c))
    return np.matrix(localXi_matrix) * np.matrix(c)


order_error_cases = np.arange(8,50)
interpolated_error = np.zeros(len(order_error_cases), dtype=np.float)

for i in range(len(order_error_cases)):
    x_axis = np.linspace(0,1,order_error_cases[i])
    no_of_test_points = 101
    x_refined = np.linspace(0,1,no_of_test_points)
    for j in range(no_of_test_points):
        interpolated_error[i] += abs((shape(x_refined[j], np.array(np.sin(4*np.pi*x_axis))))- np.sin(4*np.pi*x_refined[j]))/no_of_test_points
    print('Computing for Index = ',i)




pl.loglog(order_error_cases, interpolated_error,lw=3,label = 'Interpolated')
pl.legend()
pl.xlabel('$\mathrm{Number\;of\;points\;provided}$')
pl.ylabel('$\mathrm{L_1\;norm\;of\;error}$')
pl.title('$\mathrm{Convergence\;plot}$')
pl.loglog(order_error_cases, 15000 * (order_error_cases ** -10.999), '--', label=' $O(N^{-11})$ ')
pl.legend(loc = 'lower left').draggable()
# pl.savefig('10points.png')
pl.show()
pl.clf()

# pl.plot(x_refined, interpolated,lw=3,label = 'Interpolated')
# pl.legend()
# pl.plot(x_refined, np.array(np.sin(4*np.pi*x_refined)),'--',lw = 5,label = 'Actual')
# pl.legend()
# pl.plot(x_axis, np.array(np.sin(4*np.pi*x_axis)),'o',markersize = 25,label = 'Points Provided')
# pl.legend().draggable()
# pl.plot()
# # pl.savefig('10points.png')
# pl.show()
# pl.clf()
