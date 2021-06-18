import numpy as np
from scipy import special
import sympy as sp

## CHECK XEM HAM BOYS CO GIONG VOI TUTORIAL THU 2 KHONG
def boys(x, n):
    if x == 0:
        return 1./(2*n+1)
    else:
        return special.gammainc(n+0.5,x) * special.gamma(n+0.5) * (1./(2*x**(n+0.5)))

def gprod_1D(x1, alpha1, x2, alpha2):
    ## CHECKED TRUE
    p = alpha1 + alpha2
    q = (alpha1 * alpha2)/p
    P = (alpha1 * x1 + alpha2 * x2)/p
    Q = x1 - x2
    KAB = np.exp(-q*Q**2)
    return KAB

def ee_integral(primitive_1, primitive_2, primitive_3, primitive_4):
    aa = primitive_1.alpha
    ab = primitive_2.alpha

    c1c2 = primitive_1.coeff * primitive_2.coeff
    N = primitive_1.A * primitive_2.A
    p = primitive_1.alpha + primitive_2.alpha
    q = (primitive_1.alpha * primitive_2.alpha) / p
    Q = primitive_1.coordinates - primitive_2.coordinates
    Q2 = Q @ Q.T

    P = (primitive_1.coordinates * aa + primitive_2.coordinates * ab) / p
    A_AB = N * c1c2 * np.exp(-q*Q2)

    ########
    ac = primitive_3.alpha
    ad = primitive_4.alpha

    cc1cc2 = primitive_3.coeff * primitive_4.coeff
    NN = primitive_3.A * primitive_4.A
    pp = primitive_3.alpha + primitive_4.alpha
    qq = (primitive_3.alpha * primitive_4.alpha) / pp
    QQ = primitive_3.coordinates - primitive_4.coordinates
    QQ2 = QQ @ QQ.T

    PP = (primitive_3.coordinates * ac + primitive_4.coordinates * ad) / pp

    A_CD = NN * cc1cc2 * np.exp(-qq*QQ2)

    RPPP2 = (P-PP).T @ (P-PP)
    alpha = p*pp/(p+pp)

    primitive_ee_repulsion = A_AB * A_CD * boys(0, alpha * RPPP2) * 2 * np.pi**2.5 / (p*pp*np.sqrt(p+pp))

    return primitive_ee_repulsion

def ee_repulsion(molecule):
    n_basis = len(molecule.basis)

    V_ee = np.zeros((n_basis, n_basis, n_basis, n_basis))

    for i, orbital_1 in enumerate(molecule.basis):
        for primitive_1 in orbital_1.primitives:
            for j, orbital_2 in enumerate(molecule.basis):
                for primitive_2 in orbital_2.primitives:

                    for k, orbital_3 in enumerate(molecule.basis):
                        for primitive_3 in orbital_3.primitives:
                            for l, orbital_4 in enumerate(molecule.basis):
                                for primitive_4 in orbital_4.primitives:
                                    V_ee[i,j,k,l] += ee_integral(primitive_1,primitive_2,primitive_3,primitive_4)
    return(V_ee)