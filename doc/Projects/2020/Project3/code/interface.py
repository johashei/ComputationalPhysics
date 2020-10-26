import numpy as np
import matplotlib.pyplot as plt
import sys, os, subprocess
from mpl_toolkits import mplot3d

"""
Compile, link and run main program. Read resulting data from file and plot.
Save plot as eps with same name as the data file.
"""

# Read user input
main = sys.argv[1] # Name of the main program file
outfilename = str(input("Input name of datafile : ")) # Name of the main program's output file

# Make a directory to store data files
dir = "datafiles"
if not os.path.isdir("./"+dir):
    os.mkdir(dir)

# Compile the main program, link and run.
#''' Uncomment to not create a new data file.
subprocess.run(["c++","-c",main+".cpp"])
subprocess.run(["c++","-o",main,main+".o","physics_simulator.o","vec3.o","tests.o"],check=True)
subprocess.run(["./"+main,outfilename])
subprocess.run(["mv",outfilename,dir])  # move data file to directory'''

# Read datafile
print("Reading data from file ...",end=" ",flush=True)
infile = open(dir+"/"+outfilename,"r")
shape = eval(infile.readline())
pos = np.zeros(shape+[3])
for k,line in zip(range(int(shape[0])),infile): # k is time step number
    values = line.split(";")
    for j in range(int(shape[1])): # j is object number
        pos[k,j,:] = eval(values[j*2+1])
infile.close()
print("done")

# Plot
colors_sol = ['yellow','grey','brown','ocean blue','rust','brownish orange','light tan','aqua marine','bright blue','pale blue']
colors_esj = ['yellow','ocean blue','brownish orange']

def plot3D():
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for i in range(int(shape[1])):
        ax.plot3D(pos[:,i,0],pos[:,i,1],pos[:,i,2],'-',color="xkcd:"+colors_sol[i])

def plot2D():
    fig = plt.figure()
    ax = fig.subplots()
    for i in range(int(shape[1])):
        ax.plot(pos[:,i,0],pos[:,i,1],'-',color="xkcd:"+colors_sol[i])

    #ax.plot(pos[0,1:,0],pos[0,1:,1],'gx', pos[-1,1:,0],pos[-1,1:,1],'rx') # Start and end points
    #ax.plot(0,0,'oy') # Show sun at origin. (For when Sun is a fixed object)
    ax.set_facecolor('#010105')
    ax.set_xlabel('x [AU]',size=20)
    ax.set_ylabel('y [AU]',size=20)
    for label in(ax.get_xticklabels()+ax.get_yticklabels()):
        label.set_fontsize(20)
    plt.axis("equal")


plot2D()
plt.tight_layout()
plt.savefig('../report/'+outfilename.replace('.dat','.eps'))
plt.show()
