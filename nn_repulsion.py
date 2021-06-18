import numpy as np


def nn_repulsion(Z, atom_coordinates):
    n_atoms = len(Z)
    energy = 0

    if n_atoms == 1:
        energy = 0
    else:
        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                distance = np.linalg.norm(atom_coordinates[i] - atom_coordinates[j])
                energy = energy + Z[i] * Z[j] / distance

    return energy
