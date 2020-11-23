import matplotlib.pyplot as plt
import numpy as np
#import sys,os,subprocess

class Ising_2x2(object):
    """Analytical solution for a 2x2 lattice"""

    def __init__(self,T):
        self.T = T
        self.Z = 4*np.cosh(8/T) + 12

        self.E = (16*(np.exp(-8/self.T) - np.exp(8/self.T)))/self.Z
        self.E2 = (128*(np.exp(-8/self.T) + np.exp(8/self.T)))/self.Z
        self.abs_M =  (8*np.exp(8/self.T) + 16)/self.Z
        self.M2 = (32*np.exp(8/self.T) + 32)/self.Z
        self.C_V =  (self.E2 - self.E**2)/self.T**2
        self.chi =  (self.M2 - self.abs_M**2)/self.T

    def plot_2x2(self,ax):
        ax.plot(self.T,self.E/4,label=r'$\langle E\rangle/N$')
        ax.plot(self.T,self.abs_M/4,label=r'$\langle |M|\rangle/N$')
        ax.plot(self.T,self.E2/4)#,label=r'$\langle E^2\rangle/N$')
        ax.plot(self.T,self.M2/4)#,label=r'$\langle M^2\rangle/N$')
        ax.plot(self.T,self.C_V/4,label=r'$\langle C_V\rangle/N$')
        ax.plot(self.T,self.chi/4,label=r'$\langle \chi \rangle/N$')
        #plt.legend()

if __name__ == '__main__':
    one = Ising_2x2(1)
    print("For T = 1")
    print("E/N = {E:g} and |M|/N = {M:g}".format(E=one.E/4,M=one.abs_M/4))
    T = np.linspace(0.1,3,100)
    example = Ising_2x2(T)
    fig,ax = plt.subplots()
    example.plot_2x2(ax)
    ax.legend()
    plt.show()
