import numpy as np
import sympy as sp


class Gaussian():
    def __init__(self, alpha, coeff, coordinates=None, l1=0,l2=0,l3=0):
        self.type = "Gaussian"
        self.alpha = alpha
        self.coeff = coeff
        self.coordinates = coordinates
        self.l = (l1, l2, l2)
        # Normalization constant
        self.A = (2. * alpha / np.pi)**(0.75) # + other terms for (l1,l2,l3)

        r = sp.Symbol('r')
        self.explicit = self.A * sp.exp(-alpha*(r**2))

    def set_coordinates(self, coordinates):
        self.coordinates = np.array(coordinates)
        return self


class AtomicOrbital:
    # H_1s, H_2s
    def __init__(self, name, primitives):
        self.name = name
        self.primitives = primitives

class Atom:
    # H1, H2
    def __init__(self, name, orbitals):
        self.name = name
        self.orbitals = orbitals


class Molecule:
    def __init__(self, name, atoms):
        self.name = name
        self.atoms = atoms
        self.basis = []
        for atom in atoms:
            self.basis += atom.orbitals