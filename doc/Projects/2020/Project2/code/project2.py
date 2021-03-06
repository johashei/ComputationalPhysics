import numpy as np
import matplotlib.pyplot as plt
import sys, os, subprocess, glob

print("Compiling...")
subprocess.run(['make','main','-C','cpp_code'],check=True)
subprocess.run(['make','test','-C','cpp_code'],check=True)
subprocess.run(['mv']+glob.glob('cpp_code/*.out')+['./'],check=True)


rho_0 = float(input("input rho_0 : "))
rho_N = float(input("input rho_N : "))
N = int(input("input N : "))
algo = str(input("input algorithm : "))
M = None # Only used by Lanczos' algorithm
if (algo == "Lanczos"):
    M = int(input("input M : "))


print("Executing...")
subprocess.run(['./test.out'],check=True) # check raises error if tests fail
subprocess.run(['./main.out',str(rho_0),str(rho_N),str(N),algo,str(M)],check=True)


def analytical_solution(N):
    """Analytical solution for tridiagonal Toeplitz matrix [a,d,a]"""
    h = 1/N
    a = -1/h**2
    d = 2/h**2
    jTj = np.linspace(np.linspace(1,N-1,N-1),np.linspace(N-1,(N-1)**2,N-1),N-1)
    eigval = d + 2*a*np.cos(jTj[0,:]*np.pi*h)
    eigvec = np.sin(jTj*np.pi*h)
    eigvec = eigvec/np.sqrt(sum(eigvec[:,0]**2)) # Normalize eigen vectors
    return(eigvec,np.sort(eigval))

#[V,D] = analytical_solution(N);
#print("Analytical solution:")
#print(D[:4],"\n")
#print(V)
