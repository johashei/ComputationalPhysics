import numpy as np
import matplotlib.pyplot as plt
import subprocess

#subprocess.run(["c++","-c","physics_simulator.cpp","main.cpp"])
subprocess.run(["c++","-o","main","main.o","physics_simulator.o","vec3.o"])
subprocess.run(["./main"])

infile = open("classdata.dat","r")

shape = eval(infile.readline())

pos = np.zeros(shape+[3])

for k,line in zip(range(int(shape[0])),infile): # k is time step number
    values = line.split(";")
    for j in range(int(shape[1])): # j is object number
        pos[k,j,:] = eval(values[j*2+1])

infile.close()

plt.plot(pos[:,:,0],pos[:,:,1],'.')
plt.plot(pos[0,:,0],pos[0,:,1],'rx', pos[-1,:,0],pos[-1,:,1],'bx')
#plt.plot(rSun[:,0],rSun[:,1],'yo')
plt.axis("equal")
plt.show()
