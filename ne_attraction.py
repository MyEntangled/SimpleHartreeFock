import numpy as np
import sympy as sp
from scipy import special
from helper_fn import boys


def coulomb_attr(primitive_1, primitive_2, atom, atom_coordinates, Z):
    if primitive_1.type == "Gaussian" and primitive_2.type == 'Gaussian':
        c1c2 = primitive_1.coeff * primitive_2.coeff
        N = primitive_1.A * primitive_2.A  # normalization
        p = primitive_1.alpha + primitive_2.alpha  # auxiliary
        q = primitive_1.alpha * primitive_2.alpha / p
        Q = primitive_1.coordinates - primitive_2.coordinates
        Q2 = Q @ Q.T

        P = (primitive_1.alpha * primitive_1.coordinates) + (primitive_2.alpha * primitive_2.coordinates)
        Pp = P/p
        PG = Pp - atom_coordinates[atom]
        PG2 = PG @ PG.T

        primitive_ne_attr = -Z[atom] * N * c1c2 * np.exp(-q*Q2)  * (2.*np.pi/p) * boys(p*PG2, 0)

        return primitive_ne_attr

    else:
        pass

def ne_attraction(molecule, atom_coordinates, Z):
    n_atoms = len(Z)
    n_basis = len(molecule.basis)

    V_ne = np.zeros((n_basis, n_basis))

    for i, orbital_1 in enumerate(molecule.basis):
        for j, orbital_2 in enumerate(molecule.basis):

            for atom in range(n_atoms):
                for primitive_1 in orbital_1.primitives:
                    for primitive_2 in orbital_2.primitives:
                        V_ne[i, j] += coulomb_attr(primitive_1, primitive_2, atom, atom_coordinates, Z)
    return V_ne