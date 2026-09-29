import numpy as np

def plummer_sphere(n, M, a, G):
    m = np.full((n, 1), M / n)

    # random number generation
    x_u = np.random.uniform(0, 1, n)
    x_v = np.random.uniform(0, 1, n)
    x_w = np.random.uniform(0, 1, n)

    r = a / np.sqrt(x_u**(-2/3) - 1.0)
    x_theta = np.arccos(1.0 - 2.0 * x_v)
    x_phi = 2.0 * np.pi * x_w

    # Von Neumann rejection method
    q = np.zeros(n)
    pending_index = np.arange(n)

    remaining = n
    while remaining > 0:
        X = np.random.uniform(0, 1.0, remaining)
        Y = np.random.uniform(0, 0.1, remaining)
        
        mascara = Y < (X**2 * (1.0 - X**2)**3.5)
        valid = X[mascara]
        
        if len(valid) > 0:
            index = pending_index[:len(valid)]
            q[index] = valid
            pending_index = pending_index[len(valid):]
            remaining -= len(valid)

    # velocity magnitude expressed as a fraction of the local escape velocity
    v_mag = q * np.sqrt(2.0 * G * M / np.sqrt(r**2 + a**2))

    v_u = np.random.uniform(0, 1, n)
    v_v = np.random.uniform(0, 1, n)
    v_theta = np.arccos(1.0 - 2.0 * v_u)
    v_phi = 2.0 * np.pi * v_v

    # projection onto cartesian coordinates
    x = r * np.sin(x_theta) * np.cos(x_phi)
    y = r * np.sin(x_theta) * np.sin(x_phi)
    z = r * np.cos(x_theta)

    vx = v_mag * np.sin(v_theta) * np.cos(v_phi)
    vy = v_mag * np.sin(v_theta) * np.sin(v_phi)
    vz = v_mag * np.cos(v_theta)

    pos = np.column_stack((x, y, z))
    vel = np.column_stack((vx, vy, vz))

    return pos, vel, m