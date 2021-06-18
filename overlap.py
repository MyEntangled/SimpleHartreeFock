import numpy as np
import sympy as sp
import utility

def gprod_1D(x1, alpha1, x2, alpha2):
    ## CHECKED TRUE
    p = alpha1 + alpha2
    q = (alpha1 * alpha2)/p
    P = (alpha1 * x1 + alpha2 * x2)/p
    Q = x1 - x2
    KAB = np.exp(-q*Q**2)
    return KAB

def overlap_integral(primitive_1, primitive_2):
    if primitive_1.type == "Gaussian" and primitive_2.type == 'Gaussian':
        c1c2 = primitive_1.coeff * primitive_2.coeff
        N = primitive_1.A * primitive_2.A  # normalization
        p = primitive_1.alpha + primitive_2.alpha  # auxiliary
        q = primitive_1.alpha * primitive_2.alpha / p
        Q = primitive_1.coordinates - primitive_2.coordinates
        #
        # Q2 = Q @ Q.T
        # primitive_overlap = N * c1c2 * np.exp(-q*Q2) * (np.pi / p) ** (3/2)

        Ex = gprod_1D(primitive_1.coordinates[0], primitive_1.alpha, primitive_2.coordinates[0], primitive_2.alpha)
        Ey = gprod_1D(primitive_1.coordinates[1], primitive_1.alpha, primitive_2.coordinates[1], primitive_2.alpha)
        Ez = gprod_1D(primitive_1.coordinates[2], primitive_1.alpha, primitive_2.coordinates[2], primitive_2.alpha)
        primitive_overlap = N * c1c2 * Ex * Ey * Ez * (np.pi / p )**(3/2)
        return primitive_overlap
    else:
        pass


def overlap(molecule):
    n_basis = len(molecule.basis)
    S = np.zeros((n_basis, n_basis))

    for i, orbital_1 in enumerate(molecule.basis):
        for j, orbital_2 in enumerate(molecule.basis):

            for primitive_1 in orbital_1.primitives:
                for primitive_2 in orbital_2.primitives:
                    S[i, j] += overlap_integral(primitive_1, primitive_2)
    return S
