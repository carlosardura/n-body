import os
import h5py
import numpy as np
from plummer import plummer_sphere

if __name__ == "__main__":
    # global constants
    n = 20000
    M = 1.0
    a = 1.0
    G = 1.0

    r, v, m = plummer_sphere(n, M, a, G)

    v_cm = np.sum(m * v, axis=0) / M
    v -= v_cm

    os.makedirs('../data', exist_ok=True)
    
    with h5py.File('../data/init.h5', 'w') as f:
        f.create_dataset('positions', data=r, dtype='float64')
        f.create_dataset('velocities', data=v, dtype='float64')
        f.create_dataset('masses', data=m, dtype='float64')