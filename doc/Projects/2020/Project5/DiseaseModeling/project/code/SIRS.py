import matplotlib.pyplot as plt
import numpy as np

t = []
SIR = []
infile = open("aA.txt",'r')
for line in infile:
    words = line.split()
    #print(words[0])
    t.append(float(words[0]))
    SIR.append([float(val) for val in words[1:]])
infile.close()
t = np.array(t)
SIR = np.array(SIR)
plt.plot(t,SIR,'-')
plt.show()
