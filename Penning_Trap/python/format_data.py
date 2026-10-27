import numpy as np
from pathlib import Path

#Import and format data
def data_multi(n_particles: int, filename: Path, columns=(0, 1, 2)):
    return np.stack([
        np.loadtxt(f"{filename}{i}.txt", usecols=columns, ndmin=2)
        for i in range(n_particles)
    ])
