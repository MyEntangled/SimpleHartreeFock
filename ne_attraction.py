import numpy as np
from scipy import special
import sympy as sp

def boys(x, n):
    if x == 0:
        return 1./(2*n+1)
    else:
        return special.gammainc(n+0.5,x) * special.gamma(n+0.5) / (2*x**(n+0.5))

def gprod_1D(x1, alpha1, x2, alpha2):
    ## CHECKED TRUE
    p = alpha1 + alpha2
    q = (alpha1 * alpha2)/p
    P = (alpha1 * x1 + alpha2 * x2)/p
    Q = x1 - x2
    KAB = np.exp(-q*Q**2)
    return KAB

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


        s = N * c1c2 * np.exp(-q*Q2) * (np.pi/p)**(3/2)

        primitive_ne_attr = -Z[atom] * N * c1c2 * np.exp(-q*Q2)  * (2.*np.pi/p) * boys(p*PG2, 0)

        return primitive_ne_attr

    else:
        r = sp.Symbol('r')
        primitive_overlap = sp.integrate(primitive_1.explicit * primitive_2.explicit,
                                         (r, 0., np.inf))
        N = (sp.integrate(primitive_1.explicit * primitive_1.explicit, (r, 0., np.inf)) * sp.integrate(primitive_2.explicit * primitive_2.explicit, (r, 0., np.inf)))**(1/2)
        primitive_overlap = primitive_overlap / N

        return primitive_overlap

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
    return(V_ne)