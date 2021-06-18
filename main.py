# This is a sample Python script.

# Press ⌃R to execute it or replace it with your code.
# Press Double ⇧ to search everywhere for classes, files, tool windows, actions, and settings.

import numpy as np

from HF_driver import HF_H2_molecule

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    ## Basis set must be "STO-3G" or "6-31G"
    distance_ls = np.linspace(0.5,7.,50)
    HF_H2_molecule(distance_ls)
