# Run C++ executables, read and plot output data
import matplotlib.pyplot as plt
import numpy as np
import sys, os

try:
    sys.argv[3] == "-d"
    os.system("rm ./P1data/*") # delete datafiles from previous runs
except:
    pass

# Compile c++ code:
os.system("c++ -o Project1 project1.cpp")

# make a directory to store data:
path = "P1data"
if not os.path.isdir("./"+path):
    os.mkdir(path)

# make a figure to plot data:
fig,ax = plt.subplots(2,1)
ax[0].set_title("Numerical solution for various $n$, along with analytical solution")
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$u(x)$")
ax[1].set_title("Maximum relative error $\epsilon$ as a function of the steplength $h$")
ax[1].set_xlabel("$\log_{10}(h)$")
ax[1].set_ylabel("$\log_{10}(\epsilon)$")


# Parameters
try: # Read from command line
    algorithm = str(sys.argv[1]) # algorithm to be used by solver
    exponents = eval(sys.argv[2]) # exponent list to define size of matrix
except :
    print("Bad usage. Use %s <algorithm> <exponents>"%sys.argv[0])
    sys.exit(1)

maxerror = np.zeros(len(exponents)) # store maximum relative error
log10h = np.zeros(len(exponents)) # store log10(h)
Times = np.zeros(len(exponents)) # store timer values in seconds

for i in range(len(exponents)):
    exponent = exponents[i]
    n = int(10**exponent) # solution will be size n+2
    print("n=10^%d"%exponent)

    filename = "%s_algo_output_1e%d"%(algorithm[0],int(exponent))

    if not os.path.isfile("./"+path+"/"+filename):
        print("Data not found, calculating ...")
        # Run the compiled C++ code to get the datafile
        cmdline = "./Project1 "+filename+" "+str(n)+" "+algorithm
        cmd = cmdline
        failure = os.system(cmd)
        if failure:
            print('running Project1 failed'); sys.exit(1)

        os.system("mv "+filename+" ./"+path+"/") # place datafile in directory


    print("Reading data ...")
    infile = open("./"+path+"/"+filename,'r')

    infile.readline() # skip first line
    if int(infile.readline()) != int(n):
        print('Wrong number of steps'); sys.exit(1)

    u = np.zeros(n+2) # store numerical values
    v = np.zeros(n+2) # store analytical solution

    for j in range(n+2):
        line = infile.readline()
        values = line.split()
        u[j] = float(values[0])
        v[j] = float(values[1])
    # Read max relative error and log10(h):
    infile.readline() # skip explanation line
    errorline = infile.readline()
    values = errorline.split()
    maxerror[i] = float(values[0])
    log10h[i] = float(values[1])
    # Read timer
    infile.readline() # skip explanation line
    Times[i] = float(infile.readline())
    infile.close()

    print("Adding data to plot ...")

    if exponent > 3:
        x = np.linspace(0,1,1001) # No use plotting 10 000 points
        step = int(10**(exponent-3))
        ax[0].plot(x,u[::step],label='$n=10^%d$'%exponent)
    else:
        x = np.linspace(0,1,n+2)
        ax[0].plot(x,u,label='$n=10^%d$'%exponent)

ax[0].plot(x,v[::step],'r--',label='analytical')
ax[0].legend(loc=1)

ax[1].plot(log10h,maxerror,'*')

print("\nlog10(h) : ",log10h)
print("maxerror : ",maxerror)
print("time [s] :",Times)

plt.tight_layout()
plt.savefig("./"+path+"/P1_%d%d_%s.eps"%(exponents[0],exponents[-1],algorithm[0]))
plt.show()
