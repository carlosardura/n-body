import os
import h5py
import subprocess
import numpy as np
from plummer import plummer_sphere

if __name__ == "__main__":
    # global constants
    n = 20000
    M = 1.0
    a = 1.0
    G = 1.0

    # integration and engine parameters
    steps = 1000
    dt = 0.001
    dump_freq = 10
    theta = 0.5
    epsilon = 1e-3

    r, v, m = plummer_sphere(n, M, a, G)

    v_cm = np.sum(m * v, axis=0) / M
    v -= v_cm

    os.makedirs('../data', exist_ok=True)
    with h5py.File('../data/init.h5', 'w') as f:
        f.create_dataset('positions', data=r, dtype='float64')
        f.create_dataset('velocities', data=v, dtype='float64')
        f.create_dataset('masses', data=m, dtype='float64')

    # freezes python and hands full execution control to C++
    subprocess.run([
        "build/nbody",
        str(n),
        str(steps),
        str(dt),
        str(dump_freq),
        str(theta),
        str(epsilon)
    ], check=True) 
