import numpy as np
from scipy import special

def boys(x, n):
    if x == 0:
        return 1./(2*n+1)
    else:
        return special.gammainc(n+0.5,x) * special.gamma(n+0.5) * (1./(2*x**(n+0.5)))

def sort_eigs(eigvecs, eigvals):
    idx = eigvals.argsort()
    sorted_eigvals = eigvals[idx]
    sorted_eigvecs = eigvecs[:,idx]
    return sorted_eigvecs, sorted_eigvals