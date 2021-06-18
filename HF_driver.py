import copy

import numpy as np

from SCF import SCF
from ee_repulsion import ee_repulsion
from kinetics import kinetics
from ne_attraction import ne_attraction
from nn_repulsion import nn_repulsion
from overlap import overlap
from plot_chart import plot_chart
from utility import Gaussian, AtomicOrbital, Atom, Molecule


def HF_H2_molecule(distance_ls, basis_set=None):
    STO3G_H = [(0.3425250914E+01, 0.1543289673E+00),
               (0.6239137298E+00, 0.5353281423E+00),
               (0.1688554040E+00, 0.4446345422E+00)]
    ao_6_31G_H = [(0.1873113696E+02,0.3349460434E-01),
               (0.2825394365E+01,0.2347269535E+00),
               (0.6401216923E+00,0.8137573261E+00),
               (0.1612777588E+00,1.0000000)]

    if (basis_set is not None) and (basis_set not in ['STO-3G', '6-31G']):
        raise Exception("Basis set must be none, or either STO-3G or 6-31G!")

    if basis_set is None:
        basis_sets = ['STO-3G', '6-31G']
    else:
        basis_sets = [basis_set]

    plot_data = []

    for basis_set in basis_sets:
        energy_ls = []
        for d in distance_ls:

            if basis_set == 'STO-3G':

                ## Initialze H2 model with STO-3G
                H_pg1a = Gaussian(STO3G_H[0][0], STO3G_H[0][1])
                H_pg1b = Gaussian(STO3G_H[1][0], STO3G_H[1][1])
                H_pg1c = Gaussian(STO3G_H[2][0], STO3G_H[2][1])

                H1 = Atom('H1', [AtomicOrbital('H1_1s',
                                               [copy.deepcopy(H_pg1a).set_coordinates([0, 0, 0]),
                                                copy.deepcopy(H_pg1b).set_coordinates([0, 0, 0]),
                                                copy.deepcopy(H_pg1c).set_coordinates([0, 0, 0])])
                                 ])

                H2 = Atom('H2', [AtomicOrbital('H2_1s',
                                               [copy.deepcopy(H_pg1a).set_coordinates([d, 0, 0]),
                                                copy.deepcopy(H_pg1b).set_coordinates([d, 0, 0]),
                                                copy.deepcopy(H_pg1c).set_coordinates([d, 0, 0])])
                                 ])
                H2_molecule = Molecule("H2_mol", [H1, H2])

                Z = [1., 1.]
                N = len(Z)
                atom_coordinates = np.array([[0., 0., 0], [d, 0., 0.]])

                S = overlap(H2_molecule)
                T = kinetics(H2_molecule)
                NE_Attr = ne_attraction(H2_molecule, atom_coordinates, Z)
                EE_Repl = ee_repulsion(H2_molecule)
                NN_Repl = nn_repulsion(Z, atom_coordinates)

                H0 = T + NE_Attr
                Min_Energy = SCF(H0, EE_Repl, S, N) + NN_Repl
                print("Distance: ", d, "Energy: ", Min_Energy)
                energy_ls.append(Min_Energy)

            else:

                ## Initialize H2 model with 6-31G
                H_pg1a = Gaussian(ao_6_31G_H[0][0], ao_6_31G_H[0][1])
                H_pg1b = Gaussian(ao_6_31G_H[1][0], ao_6_31G_H[1][1])
                H_pg1c = Gaussian(ao_6_31G_H[2][0], ao_6_31G_H[2][1])
                H_pg2a = Gaussian(ao_6_31G_H[3][0], ao_6_31G_H[3][1])

                H1 = Atom('H1', [AtomicOrbital('H1_1s',
                                               [copy.deepcopy(H_pg1a).set_coordinates([0., 0., 0.]),
                                                copy.deepcopy(H_pg1b).set_coordinates([0., 0., 0.]),
                                                copy.deepcopy(H_pg1c).set_coordinates([0., 0., 0.])]),
                                 AtomicOrbital('H1_2s',
                                               [copy.deepcopy(H_pg2a).set_coordinates([0., 0., 0.])])
                                 ])

                H2 = Atom('H2', [AtomicOrbital('H2_1s',
                                               [copy.deepcopy(H_pg1a).set_coordinates([d, 0., 0.]),
                                                copy.deepcopy(H_pg1b).set_coordinates([d, 0., 0.]),
                                                copy.deepcopy(H_pg1c).set_coordinates([d, 0., 0.])]),
                                 AtomicOrbital('H2_2s',
                                               [copy.deepcopy(H_pg2a).set_coordinates([d, 0., 0.])])
                                 ])
                H2_molecule = Molecule("H2_mol", [H1, H2])

                Z = [1., 1.]
                N = len(Z)
                atom_coordinates = np.array([[0., 0., 0], [d, 0., 0.]])

                S = overlap(H2_molecule)
                T = kinetics(H2_molecule)
                NE_Attr = ne_attraction(H2_molecule, atom_coordinates, Z)
                EE_Repl = ee_repulsion(H2_molecule)
                NN_Repl = nn_repulsion(Z, atom_coordinates)

                H0 = T + NE_Attr
                Min_Energy = SCF(H0, EE_Repl, S, N) + NN_Repl

                # print("Overlap", S)
                # print("Kinetics", T)
                # print("NE_Attr", NE_Attr)
                # print("EE_Repl", EE_Repl[:,:,3,3])
                # print("NN_Repl", NN_Repl)
                print("Distance: ", d, "Energy: ", Min_Energy)
                energy_ls.append(Min_Energy)

        plot_data.append(tuple((basis_set, energy_ls)))

    plot_chart(distance_ls, plot_data)


