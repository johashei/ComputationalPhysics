import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
from scipy.optimize import curve_fit

dir = "./datafiles/"
L_list = ['40','60','80','100']

def fitfunc(T,mu,sigma,scale,const):
    # Scaled Gaussian + constant
    return(scale*norm.pdf(T,mu,sigma) + const )

figE,axE = plt.subplots()
axE.set_xlabel('T')
axE.set_ylabel(r'$\langle E \rangle/N$')

figM,axM = plt.subplots()
axM.set_xlabel('T')
axM.set_ylabel(r'$\langle |M| \rangle/N$')

figC,axC = plt.subplots()
axC.set_xlabel('T')
axC.set_ylabel(r'$\langle C_V \rangle/N$')

figX,axX = plt.subplots()
axX.set_xlabel('T')
axX.set_ylabel(r'$\langle \chi \rangle/N$')

Xtops = [] # list to fill with X top values

for L in L_list:
    infile_main = open(dir+"f"+L+".txt")
    infile_zoom = open(dir+"f"+L+"Z.txt")
    infile_cv = open(dir+"f"+L+"ZC.txt")
    infile_chi = open(dir+"f"+L+"ZX.txt")

    T = []
    expvals = []

    for infile in [infile_main,infile_zoom,infile_cv,infile_chi]:
        infile.readline()
        for line in infile:
            words = line.split()
            T.append(float(words[0]))
            expvals.append([float(val) for val in words[1:]])
        infile.close()

    expvals = np.array(expvals)
    T = np.array(T)

    colour = next(axE._get_lines.prop_cycler)['color']
    axE.plot(T[:-8],expvals[:-8,0],'.',color=colour,label=r'$%s\times%s$'%(L,L))
    axM.plot(T[:-8],expvals[:-8,1],'.',color=colour,label=r'$%s\times%s$'%(L,L))
    axC.plot(T[:-8],expvals[:-8,4],'.',color=colour,label=r'$%s\times%s$'%(L,L))
    axX.plot(T[:-8],expvals[:-8,5],'.',color=colour,label=r'$%s\times%s$'%(L,L))
    # Different marker for point with 3e6 cycles
    axE.plot(T[-8:],expvals[-8:,0],'x',color=colour)#,label=r'$%s\times%s$'%(L,L))
    axM.plot(T[-8:],expvals[-8:,1],'x',color=colour)#,label=r'$%s\times%s$'%(L,L))
    axC.plot(T[-8:],expvals[-8:,4],'x',color=colour)#,label=r'$%s\times%s$'%(L,L))
    axX.plot(T[-8:],expvals[-8:,5],'x',color=colour)#,label=r'$%s\times%s$'%(L,L))

    fitstart = 2.25#T[np.argmax(expvals[:,5])] - 0.06
    fitend = 2.36#T[np.argmax(expvals[:,5])] + 0.04
    #print(fitstart,fitend)
    Trange = np.linspace(fitstart,fitend,100)
    TtoFit = [Ti for Ti in T if fitstart<Ti<fitend]
    CtoFit = [Ci for Ci,inrange in zip(expvals[:,4],[fitstart<Ti<fitend for Ti in T]) if inrange]
    XtoFit = [Xi for Xi,inrange in zip(expvals[:,5],[fitstart<Ti<fitend for Ti in T]) if inrange]

    poptC,pcovC = curve_fit(fitfunc,TtoFit,CtoFit)
    poptX,pcovX = curve_fit(fitfunc,TtoFit,XtoFit,p0=[2.3,0.01,0.01,1])
    print(L,poptC,poptX)

    axC.plot(Trange,fitfunc(Trange,*poptC),color=colour)
    axC.axvline(poptC[0],0,1,ls='--',c=colour)
    #axX.plot(Trange,fitfunc(Trange,*poptX),color=colour)
    #axX.axvline(poptX[0],0,1,ls='--',c=colour)
    #axX.axvline(T[np.argmax(expvals[:,5])],0,1,ls='--',c=colour)
    Xtops.append(poptX[0])

Xtops = [2.320,2.300,2.295,2.290]
axE.set_prop_cycle(None)
for x in Xtops:
    colour = next(axE._get_lines.prop_cycler)['color']
    axX.axvline(x,0,1,ls='--',c=colour)

#Plot of Tc(L)
figTc,axTc = plt.subplots()
axTc.set_xlabel('L')
axTc.set_ylabel('$T_C(L)$')
L = np.array([40,60,80,100])
#TcL = Xtops
TcL = [2.28825581,2.28113579,2.27897971,2.27755436]
Lrange = np.linspace(10,150,200)
axTc.plot(L,TcL,'.')
def Tcfit(L,Tc,a):
    return(a/L + Tc)
[Tc,a],pcovL = curve_fit(Tcfit,L,TcL)
print(Tc,a)
print(2/np.log(1+np.sqrt(2)))
print((Tc - 2/np.log(1+np.sqrt(2)))/2/np.log(1+np.sqrt(2)))
axTc.plot(Lrange,a/Lrange + Tc,'-')
axTc.axhline(Tc,0,1,ls='--',c='k')

axE.legend()
axM.legend()
axC.legend()
axX.legend()
#figE.savefig('../report/energy.eps')
#figM.savefig('../report/magnetization.eps')
#figC.savefig('../report/heatcapacity.eps')
#figX.savefig('../report/susceptibility.eps')
#figTc.savefig('../report/critical.eps')
figE.tight_layout()
figM.tight_layout()
figC.tight_layout()
figX.tight_layout()
plt.show()
