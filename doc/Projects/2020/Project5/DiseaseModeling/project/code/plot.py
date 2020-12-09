import matplotlib.pyplot as plt
import numpy as np
import sys
import SIRS

infilename = sys.argv[1];
plotstep = int(sys.argv[2]);

# ODE calculation
# Setting Variables:
S_0 = 300 ; I_0 = 100 ; R_0 = 0
dt = 0.001 ; t_i = 0 ; t_f = 30
b = 2 ; beta = 0.01 ; c = 0.5

# Making plots:
plt.figure('Simulation of the spreading of a disease by the SIRS model')

plt.subplot()
S, I, R, t = SIRS.SIRSSolver(S_0,I_0,R_0, beta,b,c, dt,t_i,t_f)
SIRS.visualizer(S,I,R, t, t_unit = 'days')
plt.title("$\\beta=${beta}".format(beta=beta))

plt.subplots_adjust(left=0.12, bottom=0.10, right=0.99, top=0.95, hspace=0.40)


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

#x = np.linspace(np.log(0.001),1,100)
#plt.plot(x - np.log(0.001),np.exp(x),':r')
"""
plt.show(block=False)

Equilibration_time = float(input("Input equilibration time : "))
eqidx = np.where(t>Equilibration_time)[0][0]
print(eqidx)

MeanSIR = np.sum(SIR[eqidx:,:],0) / len(t[eqidx:])

print(MeanSIR)
"""



plt.show(block=True)
