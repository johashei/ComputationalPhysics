import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp
import sys

# NOT WORKING FIX

"""Analytical solutions for 2-dimensional Ising model, from Onsager's results
All expectation values are given in scaled units (J=1,k=1)"""

N = int(sys.argv[1])

def E(T):
    q = 2*np.sinh(2/T)/np.cosh(2/T)**2
    E = -np.cosh(2/T)/np.sinh(2/T) * \
        (1 + 2/np.pi*(2*np.tanh(2/T)**2 - 1) * sp.ellipk(q))
    return(E/N)

def M(T):
    M7N = (1 - (1 - np.tanh(1/T)**2)**4 / (16*np.tanh(1/T)**4) )**(1/8)
    return(M7N)

if __name__=="__main__":
    T = np.linspace(0.1,3,100)
    plt.plot(T,E(T),label='<E>/N')
    plt.plot(T,M(T),label='<|M|>/N')
    plt.legend()
    plt.show()
