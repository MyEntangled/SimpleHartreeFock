import numpy as np

def overlap_integral(primitive_1, primitive_2):
    if primitive_1.type == "Gaussian" and primitive_2.type == 'Gaussian':
        c1c2 = primitive_1.coeff * primitive_2.coeff
        N = primitive_1.A * primitive_2.A  # normalization
        p = primitive_1.alpha + primitive_2.alpha  # auxiliary
        q = primitive_1.alpha * primitive_2.alpha / p
        Q = primitive_1.coordinates - primitive_2.coordinates

        Q2 = Q @ Q.T
        primitive_overlap = N * c1c2 * np.exp(-q*Q2) * (np.pi / p) ** (3/2)

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
