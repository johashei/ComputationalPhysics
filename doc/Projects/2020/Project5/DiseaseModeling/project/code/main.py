import matplotlib.pyplot as plt
import numpy as np
import sys
import subprocess
import SIRS

# Basic SIRS parameters:
a0 = float(sys.argv[1])
b = float(sys.argv[2])
c = float(sys.argv[3])
S = int(sys.argv[4])
I = int(sys.argv[5])
R = int(sys.argv[6])
Tfinal = float(sys.argv[7])
# Vital dynamics
d = float(sys.argv[8])
dI = float(sys.argv[9])
e = float(sys.argv[10])
# Seasonal variation
A = float(sys.argv[11])
omega = float(sys.argv[12])
# Vaccination
f0 = float(sys.argv[13])



# Solve with SIRS.py
def f(t):
    #return f0*(20<t<40) # Campaign
    return f0*(t-20)*(20<t) # Increasing rate

solver = SIRS.SIRS(S,I,R)
solver.Basic_SIRS(a0,b,c)
solver.Vital_dynamics(e,d,dI)
solver.Seasonal_variation(A,omega)
solver.Vaccination(f)
SIR_ode,t_ode = solver.Solve_ODE(Tfinal,0.01)


# Compile and run main.cpp with the same argumentss
subprocess.run(['c++','-o','MCmain','main.cpp','SIRS.o','-std=c++11'])
cmd = ['./MCmain']+[str(x) for x in [a0,b,c,S,I,R,Tfinal,d,dI,e,A,omega,f0]]
subprocess.run(cmd)
# Read the results file
t_mc = []
SIR_mc = []
infile = open("SIRSmc.tsv",'r')
for line in infile:
    words = line.split()
    #print(words[0])
    t_mc.append(float(words[0]))
    SIR_mc.append([float(val) for val in words[1:]])
infile.close()
t_mc = np.array(t_mc)
SIR_mc = np.array(SIR_mc)


def MakePlots():
    colours = ['xkcd:blue','xkcd:vomit','xkcd:hospital green']
    labels = ['Susceptible','Infected','Recovered']
    fig,ax = plt.subplots(figsize=[4.5,4])
    for i in range(3):
        ax.plot(t_ode,SIR_ode[:,i],c=colours[i],ls='-',label=labels[i])
        ax.step(t_mc,SIR_mc[:,i],c=colours[i],ls='-')
    #ax.plot(t_ode,np.sum(SIR_ode,axis=1),'m-')
    #ax.step(t_mc,np.sum(SIR_mc,axis=1),'m-')
    ax.set_xlabel('time (arbitrary units)')
    ax.set_ylabel('population')

    ax.legend()
    fig.tight_layout()
    return(fig,ax)

fig,ax = MakePlots()
fig.show()


# Get equilibration time and calculate mean and standard deviation
try:
    Equilibration_time = float(input("Input equilibration time : "))
except:
    plt.close(fig)
    fig,ax = MakePlots()
    plt.show()

plt.close(fig)
eqidx = np.where(t_mc>Equilibration_time)[0][0]

MeanSIR = np.mean(SIR_mc[eqidx:,:],axis=0)
sigmaSIR = np.std(SIR_mc[eqidx:,:],axis=0)

#print("After equilibration time", Equilibration_time)
print(r"  &    Monte Carlo     &  ODE  \\")
for i in range(3):
    print(r"{C} & {mu:6.1f} $\pm$ {sigma:5.1f} & {ode:6.1f} \\"\
    .format(mu=MeanSIR[i],sigma=sigmaSIR[i],C='SIR'[i],ode=SIR_ode[-1,i]))

fig,ax = MakePlots()

for mu,sigma in zip(MeanSIR,sigmaSIR):
    ax.plot([t_mc[eqidx],t_mc[-1]],[mu,mu],c='r',ls='--')
    ax.plot([t_mc[eqidx],t_mc[-1]],[mu+sigma,mu+sigma],c='r',ls=':')
    ax.plot([t_mc[eqidx],t_mc[-1]],[mu-sigma,mu-sigma],c='r',ls=':')

plt.show()
