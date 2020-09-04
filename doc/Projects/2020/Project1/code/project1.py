# Run C++ executables, read and plot output data
import matplotlib.pyplot as plt
import numpy as np
import sys, os

os.system("make reset") # delete datafiles from previous runs

# Compile c++ code:
os.system("c++ -o Project1b Project1b.cpp")

# make a directory to store data:
path = "P1data"
if not os.path.isdir("./"+path):
    os.mkdir(path)

# Parameters
exponents = np.linspace(1,5,5) # exponents define size of matrix
algorithm = 'specific'

maxerror = np.zeros(len(exponents)) # store maximum relative error
log10h = np.zeros(len(exponents)) # store log10(h)
Times = np.zeros(len(exponents)) # store timer values in s

for i in range(len(exponents)):
    exponent = exponents[i]
    n = int(10**exponent) # solution will be size n+2
    print("n=10^%d"%exponent)

    filename = "%s_algo_output_1e%d"%(algorithm[0],int(exponent))

    if not os.path.isfile("./"+path+"/"+filename):
        print("Data not found, calculating ...")
        # Run the compiled C++ code to get the datafile
        cmdline = "./Project1b "+filename+" "+str(n)+" "+algorithm
        cmd = cmdline
        failure = os.system(cmd)
        if failure:
            print('running project1 failed'); sys.exit(1)

        os.system("mv "+filename+" ./"+path+"/") # place datafile in directory
        print("Done!")


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

    print("Adding data to plot")

    x = np.linspace(0,1,n+2)
    plt.plot(x,u,label='$n=10^%d$'%exponent)

plt.plot(x,v,'r--',label='analytical')
plt.legend(loc=1)

print("\nlog10(h) : ",log10h)
print("maxerror : ",maxerror)
print("time [s] :",Times)

plt.show()
