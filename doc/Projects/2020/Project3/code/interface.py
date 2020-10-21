import numpy as np
import matplotlib.pyplot as plt
import subprocess
from mpl_toolkits import mplot3d

#subprocess.run(["c++","-c","physics_simulator.cpp","main.cpp"])
subprocess.run(["c++","-o","main","main.o","physics_simulator.o","vec3.o"],check=True)
subprocess.run(["./main"])

print("Reading data from file ...",end=" ",flush=True)
infile = open("classdata.dat","r")

shape = eval(infile.readline())

pos = np.zeros(shape+[3])

for k,line in zip(range(int(shape[0])),infile): # k is time step number
    values = line.split(";")
    for j in range(int(shape[1])): # j is object number
        pos[k,j,:] = eval(values[j*2+1])

infile.close()
print("done")

colors = ['yellow','grey','brown','ocean blue','rust','brownish orange','light tan','aqua marine','bright blue','pale blue']

fig = plt.figure()
ax = plt.axes(projection='3d')
for i in range(int(shape[1])):
    ax.plot3D(pos[:,i,0],pos[:,i,1],pos[:,i,2],'-',color="xkcd:"+colors[i])

#ax.plot(pos[0,:,0],pos[0,:,1],'rx', pos[-1,:,0],pos[-1,:,1],'bx')
#ax.set_facecolor('xkcd:dark grey')
#plt.plot(rSun[:,0],rSun[:,1],'yo')
#plt.axis("equal")
plt.show()
