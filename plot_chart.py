import matplotlib.pyplot as plt
plt.style.use('ggplot')

def plot_chart(x,y,basis_set):
    print(len(x), len(y))
    plt.plot(x,y, label=basis_set)
    plt.xlabel("Distance (a.u.)")
    plt.ylabel("Energy (Hartree)")
    plt.title("H2 bond dissociation curve")
    plt.legend()
    plt.show()