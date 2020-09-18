import numpy as np
import matplotlib.pyplot as plt
import sys, os, subprocess, glob

N = int(input("input N : "))
filename = "Jacobirot_N_{N:d}.txt".format(N=N)

path_to_data = "./data"
if not os.path.isdir(path_to_data):
    os.mkdir(path_to_data)

# Compile C++ code:
print("Compiling...")
all_cpp_codes = glob.glob("./cpp_code/*.cpp")
subprocess.run(['c++','-o','main.out']+all_cpp_codes+['-larmadillo'], check=True)

print("Executing...")
subprocess.run(['./main.out',str(N),filename])

subprocess.run(['mv',filename,path_to_data])

def analytical_solution(N): # Should be equal to numerical to a factor
    h = 1/N
    a = -1/h**2
    d = 2/h**2
    jTj = np.linspace(np.linspace(1,N-1,N-1),np.linspace(N-1,(N-1)**2,N-1),N-1)
    eigval = d + 2*a*np.cos(jTj[0,:]*np.pi*h)
    eigvec = np.sin(jTj*np.pi*h)
    return(eigvec,eigval)

[V,D] = analytical_solution(N);
print(D)
print(V)
