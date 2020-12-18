import matplotlib.pyplot as plt
import numpy as np
import sys
import SIRS

infilename = sys.argv[1];
plotstep = int(sys.argv[2]);


t = []
SIR = []
infile = open(infilename,'r')
for line in infile:
    words = line.split()
    #print(words[0])
    t.append(float(words[0]))
    SIR.append([float(val) for val in words[1:]])
infile.close()
t = np.array(t)
SIR = np.array(SIR)
#N = np.sum(SIR[0,:])
plt.step(t[::plotstep],SIR[::plotstep,:],'-')
#plt.step(t[::plotstep],N-np.sum(SIR,1),'-')



plt.show()
