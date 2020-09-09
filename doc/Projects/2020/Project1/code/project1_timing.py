import numpy as np
import sys, os, subprocess

# Program for timing algorithms in project1.cpp. Does not store algo results.

# Compile C++ code:
os.system("c++ -o Project1 project1.cpp lib.cpp")

# Read parameters from command line:
try:
    algorithm = str(sys.argv[1]) # algorithm to be used by solver
    n = int(sys.argv[2]) # exponent list to define size of matrix
    sample_size = int(sys.argv[3]) # sample size for timing algorithm
except :
    print("Bad usage. Use %s <algorithm> <number of steps> <sample size>"%sys.argv[0])
    sys.exit(1)

# Initialise array:
Times = np.zeros(sample_size)

for i in range(sample_size):
    # Run C++ code and read output as string:
    cmd = "./Project1 None "+str(n)+" "+algorithm
    output = subprocess.run([cmd], stdout=subprocess.PIPE, shell=True).stdout.decode('utf-8')
    # Find time value and convert to float
    Times[i] = float(output[output.index('{')+1 : output.index("}")])

avg_time = np.mean(Times)

print("For %s algo with %d steps, mean runtime over %d runs:"%(algorithm,n,sample_size))
print(avg_time)
