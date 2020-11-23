import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import sys,os,subprocess

dir = 'datafiles/'
names = ['d1a.txt','d24r.txt']

fig,ax = plt.subplots(2,1)

start = int(5e5)

for i in range(2):
    name = names[i]
    # Read file
    infile = open(dir+name)
    E = []
    words = infile.readline().split()
    L = int(words[0])
    MCS = int(words[1])
    for line in infile:
        words = line.split()
        E.append(float(words[0]))
    infile.close()

    # Plot histogram
    E = E[int(start):]
    MCS = MCS-int(start)
    N = L**2

    exp_E = np.sum(np.array(E))/MCS
    exp_E2 = np.sum(np.array(E)**2)/MCS

    Nbins = int((max(E)-min(E))/4) # bin width of 4
    ax[i].hist(E,Nbins,density=True)
    ax[i].set_xlim((min(E),max(E)))

    sigmaE2 = exp_E2 - exp_E**2
    sigmaE = np.sqrt(sigmaE2)
    x = np.linspace(min(E),max(E),100)
    print(exp_E,sigmaE)
    ax[i].plot(x,stats.norm.pdf(x,exp_E,sigmaE))
    ax[i].set_xlabel("E")
    ax[i].set_ylabel("P(E)")
    ax[i].ticklabel_format(axis='y',style='sci',scilimits=(0,0),useMathText=True)

plt.tight_layout()
plt.show()
