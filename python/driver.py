import h5py
import subprocess
import numpy as np
from pathlib import Path
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

    root = Path(__file__).resolve().parent.parent
    
    init_file_path = root / 'data' / 'init.h5'
    output_dir_path = root / 'data' / 'output'
    nbody_executable = root / 'build' / 'nbody'

    r, v, m = plummer_sphere(n, M, a, G)

    v_cm = np.sum(m * v, axis=0) / M
    v -= v_cm

    output_dir_path.mkdir(parents=True, exist_ok=True)
    with h5py.File(init_file_path, 'w') as f:
        f.create_dataset('positions', data=r, dtype='float64')
        f.create_dataset('velocities', data=v, dtype='float64')
        f.create_dataset('masses', data=m, dtype='float64')

    # freezes python and hands full execution control to C++
    subprocess.run([
        str(nbody_executable),
        str(init_file_path),
        str(output_dir_path),
        str(n),
        str(G),
        str(steps),
        str(dt),
        str(dump_freq),
        str(theta),
        str(epsilon)
    ], check=True)