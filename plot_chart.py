import matplotlib.pyplot as plt
import numpy as np
plt.style.use('ggplot')

def plot_chart(x,y):
    basis_sets, energy_sets = zip(*y)
    for idx in range(len(basis_sets)):
        plt.plot(x,energy_sets[idx], label=basis_sets[idx])

    plt.xlabel("Distance (a.u.)")
    plt.ylabel("Energy (Hartree)")
    plt.xticks(np.arange(x[0], x[-1], step=0.5))
    plt.title("H2 bond dissociation curve")
    plt.legend()
    plt.savefig('combined_result.png')
    plt.show()