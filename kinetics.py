import numpy as np

def kinetics_integral(primitive_1, primitive_2):
    if primitive_1.type == "Gaussian" and primitive_2.type == 'Gaussian':
        c1c2 = primitive_1.coeff * primitive_2.coeff
        N = primitive_1.A * primitive_2.A  # normalization
        p = primitive_1.alpha + primitive_2.alpha  # auxiliary
        q = primitive_1.alpha * primitive_2.alpha / p
        Q = primitive_1.coordinates - primitive_2.coordinates
        Q2 = Q @ Q.T

        P = (primitive_1.alpha * primitive_1.coordinates) + (primitive_2.alpha * primitive_2.coordinates)
        Pp = P/p
        PG = Pp - primitive_2.coordinates
        PGx2 = PG[0] * PG[0]
        PGy2 = PG[1] * PG[1]
        PGz2 = PG[2] * PG[2]

        s = N * c1c2 * np.exp(-q*Q2) * (np.pi/p)**(3/2)
        primitive_kinetics = 3*primitive_2.alpha * s \
                            -2*primitive_2.alpha * primitive_2.alpha * s * (PGx2 + 0.5/p) \
                            -2*primitive_2.alpha * primitive_2.alpha * s * (PGy2 + 0.5/p) \
                            -2*primitive_2.alpha * primitive_2.alpha * s * (PGz2 + 0.5/p)

        return primitive_kinetics

    else:
        pass


def kinetics(molecule):
    n_basis = len(molecule.basis)
    T = np.zeros((n_basis, n_basis))

    for i, orbital_1 in enumerate(molecule.basis):
        for j, orbital_2 in enumerate(molecule.basis):
            T[i, j] = 0

            for primitive_1 in orbital_1.primitives:
                for primitive_2 in orbital_2.primitives:
                    T[i, j] += kinetics_integral(primitive_1, primitive_2)

    return T
