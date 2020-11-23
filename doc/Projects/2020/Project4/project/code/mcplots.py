import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import sys,os,subprocess

dir = 'datafiles/'
names = ['d1a.txt','d1r.txt','d24a.txt','d24r.txt']

fig,ax = plt.subplots(3,1,sharex=True)
#figE,ax[0] = plt.subplots()
#figM,ax[1] = plt.subplots()
#figNa,ax[2] = plt.subplots()

for name in names:
    # Read file
    infile = open(dir+name)
    E = []
    M = []
    N_accepted_configs = []
    words = infile.readline().split()
    L = int(words[0])
    MCS = int(words[1])
    T = float(words[2])
    spin0 = words[3]
    for line in infile:
        words = line.split()
        E.append(float(words[0]))
        M.append(float(words[1]))
        N_accepted_configs.append(int(words[2]))
    infile.close()

    # Calculate values and plot
    MCS = np.linspace(1,MCS,MCS)
    N = L**2

    exp_E = np.cumsum(np.array(E))/MCS
    exp_E2 = np.cumsum(np.array(E)**2)/MCS
    exp_M = np.cumsum(abs(np.array(M)))/MCS
    exp_M2 = np.cumsum(np.array(M)**2)/MCS
    N_total = np.cumsum(N_accepted_configs)/MCS

    X = []
    YE,YM,YN = [],[],[]
    idx0 = 0
    # Plot at most 100 points for each power of 10
    for i in range(1,int(np.ceil(np.log10(MCS[-1])))+1):
        idxf = 10**i
        step = max([int(idxf/30),1])
        X += list(MCS[idx0:idxf:step])
        YE += list(exp_E[idx0:idxf:step]/N)
        YM += list(exp_M[idx0:idxf:step]/N)
        YN += list(np.array(N_accepted_configs)[idx0:idxf:step]/N)
        idx0 = idxf

    ax[0].semilogx(X,YE,'.-',label=name[1:-4])
    ax[1].semilogx(X,YM,'.-',label=name[1:-4])
    ax[2].semilogx(X,YN,'.-',label=name[1:-4])


#ax[0].set_xlabel("Monte Carlo cycles")
ax[0].set_ylabel(r"$\langle E\rangle/N$")
#ax[0].legend()

#ax[1].set_xlabel("Monte Carlo cycles")
ax[1].set_ylabel(r"$\langle |M|\rangle/N$")
#ax[1].legend()

ax[2].set_xlabel("Monte Carlo cycles")
ax[2].set_ylabel("A/N")
ax[0].legend(loc='upper right',ncol=4)

plt.show()
