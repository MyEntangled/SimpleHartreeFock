# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np
import copy

from utility import Gaussian, AtomicOrbital, Atom, Molecule
from overlap import overlap
from kinetics import kinetics
from ne_attraction import ne_attraction
from ee_repulsion import ee_repulsion
from nn_repulsion import nn_repulsion
from SCF import SCF
from plot_chart import plot_chart
from HF_driver import HF_H2_molecule

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    ## Basis set must be "STO-3G" or "6-31G"
    distance_ls = np.linspace(0.6,7,50)
    HF_H2_molecule("6-31G", distance_ls)

    # # STO3G_H = [(0.3425250914E+01, 0.1543289673E+00),
    # #            (0.6239137298E+00, 0.5353281423E+00),
    # #            (0.1688554040E+00, 0.4446345422E+00)]
    # # H_pg1a = Gaussian(STO3G_H[0][0], STO3G_H[0][1])
    # # H_pg1b = Gaussian(STO3G_H[1][0], STO3G_H[1][1])
    # # H_pg1c = Gaussian(STO3G_H[2][0], STO3G_H[2][1])
    # ao_4_31G_H = [(0.1873113696E+02,0.3349460434E-01),
    #            (0.2825394365E+01,0.2347269535E+00),
    #            (0.6401216923E+00,0.8137573261E+00),
    #            (0.1612777588E+00,1.0000000)]
    # H_pg1a = Gaussian(ao_4_31G_H[0][0], ao_4_31G_H[0][1])
    # H_pg1b = Gaussian(ao_4_31G_H[1][0], ao_4_31G_H[1][1])
    # H_pg1c = Gaussian(ao_4_31G_H[2][0], ao_4_31G_H[2][1])
    # H_pg2a = Gaussian(ao_4_31G_H[3][0], ao_4_31G_H[3][1])
    #
    # Energy_ls = []
    # dist = np.linspace(0.6, 7., 30)
    # for d in dist:
    #     #print(d)
    #     # H1 = Atom('H1', [AtomicOrbital('H1_1s',
    #     #                                [copy.deepcopy(H_pg1a).set_coordinates([0, 0, 0]),
    #     #                                 copy.deepcopy(H_pg1b).set_coordinates([0, 0, 0]),
    #     #                                 copy.deepcopy(H_pg1c).set_coordinates([0, 0, 0])])
    #     #                  ])
    #     #
    #     # H2 = Atom('H2', [AtomicOrbital('H2_1s',
    #     #                                [copy.deepcopy(H_pg1a).set_coordinates([d, 0, 0]),
    #     #                                 copy.deepcopy(H_pg1b).set_coordinates([d, 0, 0]),
    #     #                                 copy.deepcopy(H_pg1c).set_coordinates([d, 0, 0])])
    #     #                  ])
    #     H1 = Atom('H1', [AtomicOrbital('H1_1s',
    #                                    [copy.deepcopy(H_pg1a).set_coordinates([0., 0., -d/2]),
    #                                     copy.deepcopy(H_pg1b).set_coordinates([0., 0., -d/2]),
    #                                     copy.deepcopy(H_pg1c).set_coordinates([0., 0., -d/2])]),
    #                      AtomicOrbital('H1_2s',
    #                                    [copy.deepcopy(H_pg2a).set_coordinates([0., 0., -d/2])])
    #                      ])
    #
    #     H2 = Atom('H2', [AtomicOrbital('H2_1s',
    #                                    [copy.deepcopy(H_pg1a).set_coordinates([0., 0., d/2]),
    #                                     copy.deepcopy(H_pg1b).set_coordinates([0., 0., d/2]),
    #                                     copy.deepcopy(H_pg1c).set_coordinates([0., 0., d/2])]),
    #                      AtomicOrbital('H2_2s',
    #                                    [copy.deepcopy(H_pg2a).set_coordinates([0., 0., d/2])])
    #                      ])
    #     Z = [1., 1.]
    #     atom_coordinates = np.array([[0., 0., -d/2], [0, 0., d/2]])
    #     H2_molecule = Molecule("H2_mol", [H1, H2])
    #
    #     S = overlap(H2_molecule)
    #     T = kinetics(H2_molecule)
    #     NE_Attr = ne_attraction(H2_molecule, atom_coordinates, Z)
    #     EE_Repl = ee_repulsion(H2_molecule)
    #
    #     ## N = TOTAL NUMBER OF ELECTRONS IN THE SYSTEM
    #     ## H0: SINGLE ELECTRON HAMILTONIAN
    #
    #     N = len(Z)
    #     H0 = T + NE_Attr
    #     Min_Energy = SCF(H0, EE_Repl, S, N)
    #     Min_Energy += nn_repulsion(Z,atom_coordinates)
    #
    #     Energy_ls.append(Min_Energy)
    #
    # plot_chart(dist, Energy_ls, "6-31G")
    #
    #
    # print("STO-3G for H2")
    # print("Overlap", S)
    # print("Kinetics", T)
    # print("NE Attraction", NE_Attr)
    # print("EE Repulsion", EE_Repl[:,:,-1,-1])
    # print("MIN ENERGY", Min_Energy)
    #
    #
    # bs_6_31G = [(0.1873113696E+02, 0.3349460434E-01),
    #             (0.2825394365E+01, 0.2347269535E+00),
    #             (0.6401216923E+00, 0.8137573261E+00),
    #             (0.1612777588E+00, 1.0000000)]
    #
    # H_pg1a = Gaussian(bs_6_31G[0][0], bs_6_31G[0][1])
    # H_pg1b = Gaussian(bs_6_31G[1][0], bs_6_31G[1][1])
    # H_pg1c = Gaussian(bs_6_31G[2][0], bs_6_31G[2][1])
    # H_pg2a = Gaussian(bs_6_31G[3][0], bs_6_31G[3][1])
    #
    # H1 = Atom('H1', [AtomicOrbital('H1_1s',
    #                                [copy.deepcopy(H_pg1a).set_coordinates([0, 0, 0]),
    #                                 copy.deepcopy(H_pg1b).set_coordinates([0, 0, 0]),
    #                                 copy.deepcopy(H_pg1c).set_coordinates([0, 0, 0])]),
    #                  AtomicOrbital('H1_2s',
    #                                [copy.deepcopy(H_pg2a).set_coordinates([0, 0, 0])])
    #                  ])
    #
    # H2 = Atom('H2', [AtomicOrbital('H2_1s',
    #                                [copy.deepcopy(H_pg1a).set_coordinates([d, 0, 0]),
    #                                 copy.deepcopy(H_pg1b).set_coordinates([d, 0, 0]),
    #                                 copy.deepcopy(H_pg1c).set_coordinates([d, 0, 0])]),
    #                  AtomicOrbital('H2_2s',
    #                                [copy.deepcopy(H_pg2a).set_coordinates([d, 0, 0])])
    #                  ])
    # Z = [1., 1.]
    # atom_coordinates = np.array([[0.,0.,0.], [d, 0, 0]])
    # H2_molecule = Molecule("H2_mol", [H1, H2])
    # print("6-31G for H2")
    # print("Overlap", overlap(H2_molecule))
    # print("Kinetics", kinetics(H2_molecule))
    # print("NE Attraction", ne_attraction(H2_molecule, atom_coordinates, Z))
    # print("EE Repulsion", ee_repulsion(H2_molecule).shape)
