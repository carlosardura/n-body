import numpy as np

def plummer_sphere(n, M, a):
    m = np.full((n, 1), M / n)

    # random number generation
    x_u = np.random.uniform(0, 1, n)
    x_v = np.random.uniform(0, 1, n)
    x_w = np.random.uniform(0, 1, n)

    r = a / np.sqrt(x_u**(-2/3) - 1.0)
    x_theta = np.arccos(1.0 - 2.0 * x_v)
    x_phi = 2.0 * np.pi * x_w

    # projection onto cartesian coordinates
    x = r * np.sin(x_theta) * np.cos(x_phi)
    y = r * np.sin(x_theta) * np.sin(x_phi)
    z = r * np.cos(x_theta)

    pos = np.column_stack((x, y, z))

    return pos, m