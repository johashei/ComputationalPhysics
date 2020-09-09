import numpy as np
import sys, os, subprocess

# Program for timing algorithms in project1.cpp. Does not store algo results.

# Compile C++ code:
os.system("c++ -o Project1 project1.cpp lib.cpp")

# Read parameters from command line:
try:
    algorithm = str(sys.argv[1]) # algorithm to be used by solver
    exponents = eval(sys.argv[2]) # exponent list to define size of matrix
    sample_size = int(sys.argv[3]) # sample size for timing algorithm
except :
    print("Bad usage. Use %s <algorithm> <exponents> <sample size>"%sys.argv[0])
    sys.exit(1)

# Initialise array:
Times = np.zeros([len(exponents),sample_size])

for i in range(len(exponents)):
    exponent = exponents[i]
    n = int(10**exponent)
    for j in range(sample_size):
        # Run C++ code and read output as string:
        cmd = "./Project1 None "+str(n)+" "+algorithm
        output = subprocess.run([cmd], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
        # Find time value and convert to float
        Times[i][j] = float(output[output.index('{')+1 : output.index("}")])

# Compute mean and standard deviation along sample axis
avg_times = np.mean(Times,1)
err_times = np.std(Times,1)/np.sqrt(sample_size)


# Print result in LaTeX tabular format
print(" Whith %s algorithm, over %d runs:"%(algorithm,sample_size))
print(r" n  &  mean : error in the mean \\")
for i in range(len(exponents)):
    print(r"1e{e:1d} & {m:10.5f} : {s:8.5f} \\".format(e=int(exponents[i]),m=avg_times[i],s=err_times[i]))
