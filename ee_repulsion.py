import numpy as np
from helper_fn import boys


def one_pair_integral(primitive_a, primitive_b):
    c1c2 = primitive_a.coeff * primitive_b.coeff
    N = primitive_a.A * primitive_b.A
    p = primitive_a.alpha + primitive_b.alpha
    q = (primitive_a.alpha * primitive_b.alpha) / p
    Q = primitive_a.coordinates - primitive_b.coordinates
    Q2 = Q @ Q.T

    P = (primitive_a.coordinates * primitive_a.alpha + primitive_b.coordinates * primitive_b.alpha) / p
    A_AB = N * c1c2 * np.exp(-q * Q2)

    return p,P,A_AB


def two_pair_integral(pair_1_res, pair_2_res):
    p,P,A_AB = pair_1_res
    pp,PP,A_CD = pair_2_res

    RPPP2 = (P-PP) @ (P-PP).T
    alpha = p*pp/(p+pp)
    primitive_ee_repulsion = A_AB * A_CD * boys(alpha * RPPP2, 0) * 2 * np.pi ** 2.5 / (p * pp * np.sqrt(p + pp))
    return primitive_ee_repulsion

def ee_repulsion(molecule):
    n_basis = len(molecule.basis)

    V_ee = np.zeros((n_basis, n_basis, n_basis, n_basis))

    for i, orbital_1 in enumerate(molecule.basis):
        for primitive_1 in orbital_1.primitives:
            for j, orbital_2 in enumerate(molecule.basis):
                for primitive_2 in orbital_2.primitives:

                    pair_1_res = one_pair_integral(primitive_1, primitive_2)

                    for k, orbital_3 in enumerate(molecule.basis):
                        for primitive_3 in orbital_3.primitives:
                            for l, orbital_4 in enumerate(molecule.basis):
                                for primitive_4 in orbital_4.primitives:

                                    pair_2_res = one_pair_integral(primitive_3, primitive_4)
                                    V_ee[i,j,k,l] += two_pair_integral(pair_1_res, pair_2_res)
    return V_ee