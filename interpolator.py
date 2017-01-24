import numpy as np
import numpy.linalg as la


"""------------------------Basic-Bilinear-Interpolation-Function(Scipy libraries can also be used)------------------ """


def bilinear_interpolate(x, y, x_grid, y_grid, F):
    F_function = F

    n = len(F_function[0, :]) - 3  # number of zones

    x_zone = int(n * (x - x_grid[0]))  # indexing from zero itself
    y_zone = int(n * (y - y_grid[0]))

    b = np.matrix([[0], [0], [0], [0]])

    A = np.matrix( \
        [[1, x_grid[x_zone], y_grid[y_zone], x_grid[x_zone] * y_grid[y_zone]], \
         [1, x_grid[x_zone], y_grid[y_zone + 1], x_grid[x_zone] * y_grid[y_zone + 1]], \
         [1, x_grid[x_zone + 1], y_grid[y_zone], x_grid[x_zone + 1] * y_grid[y_zone]], \
         [1, x_grid[x_zone + 1], y_grid[y_zone + 1], x_grid[x_zone + 1] * y_grid[y_zone + 1]] \
         ]
    )

    point_to_calculated_for = np.matrix([[1], [x], [y], [x * y]])

    b = (la.inv(A)).transpose() * point_to_calculated_for

    Q11 = F_function[y_zone, x_zone]
    Q21 = F_function[y_zone, x_zone + 1]
    Q12 = F_function[y_zone + 1, x_zone]
    Q22 = F_function[y_zone + 1, x_zone + 1]

    Q = np.matrix([[Q11], [Q12], [Q21], [Q22]])

    F_interpolated = b.transpose() * Q

    return F_interpolated


bilinear_interpolate = np.vectorize(bilinear_interpolate, excluded=(['x_grid', 'y_grid', 'F']))

"""-------------------------------------------------------END--------------------------------------------------------"""
