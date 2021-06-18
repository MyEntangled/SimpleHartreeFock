import numpy as np
from scipy.linalg import fractional_matrix_power
from helper_fn import sort_eigs

def fock_energy(D, H0, F):
    dim = D.shape[0]
    E = 0

    for n in range(dim):
        for m in range(dim):
            E += D[n][m]*(H0[n][m] + F[n][m])

    return E

def Build_Coulomb_Exchange(D, EE_Repl):
    dim = D.shape[0]
    G = np.zeros((dim, dim))
    for i in range(dim):
        for j in range(dim):
            for k in range(dim):
                for l in range(dim):
                    ## Swaping basis set 2 and 3 to get exchange integral
                    G[i,j] = G[i,j] + D[k,l] * (2*EE_Repl[i,j,k,l] - EE_Repl[i,k,j,l])

    return G

def Build_Density(eigvecs, N):
    dim = eigvecs.shape[0]
    D = np.zeros((dim,dim))

    for n in range(dim):
        for m in range(dim):
            for i in range(int(N/2)):
                D[n,m] = D[n,m] + eigvecs[n,i] * eigvecs[m,i]

    return D

def SCF(H0, EE_Repl, S, N_elec):
    maxcycles = 30
    converged = 0
    ncycle = 0
    E = np.zeros(maxcycles)

    eigvals, V = np.linalg.eigh(S)
    D = np.diag(eigvals)

    #X = V @ fractional_matrix_power(D, -0.5) @ V.T
    X = V @ fractional_matrix_power(D, -0.5)

    F = H0
    Fprime = X.T @ F @ X

    eps, Cprime = np.linalg.eig(Fprime)
    C = X @ Cprime

    C, eps = sort_eigs(C, eps)
    D = Build_Density(C,N_elec)

    while ncycle < maxcycles-1 and converged != 1:
        ncycle += 1
        G = Build_Coulomb_Exchange(D, EE_Repl)
        F = H0 + G
        Fprime = X.T @ F @ X

        eps, Cprime = np.linalg.eig(Fprime)
        C = X @ Cprime

        C, eps = sort_eigs(C, eps)
        D = Build_Density(C, N_elec)

        E[ncycle] = fock_energy(D, H0, F)

        if (ncycle > 1) and (abs(E[ncycle] - E[ncycle-1]) < 10e-9):
            converged = 1
    MinE = E[ncycle]
    return MinE