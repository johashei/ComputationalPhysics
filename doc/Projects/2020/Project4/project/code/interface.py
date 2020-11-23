import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import sys,os,subprocess
import Ising_2x2

try:
    outfilename = sys.argv[1]
except:
    run_sim = str(input("run simulation ? [y/n] : "))

    if run_sim == "y":
        print('\033[1m'+"Input simulation parameters"+'\033[0m')
        L = str(input("Side length of the lattice : "))
        print("Monte Carlo cycles")
        MCSi = str(input("\tlog10(Start) : "))
        MCSf = str(input("\tlog10(Stop) : "))
        dMCS = str(input("\tlog10(Step) : "))
        print("Temperature")
        Ti = str(input("\tStart : "))
        Tf = str(input("\tStop : "))
        dT = str(input("\tStep : "))
        spini = str(input("Initial spin state : "))
        outfilename = str(input("Name of file to generate : "))

        subprocess.run(['./main_omp',L,MCSi,MCSf,dMCS,Ti,Tf,dT,outfilename,spini],check=True)
    elif run_sim == "n":
        outfilename = str(input("Name of datafile to plot : "))


def plot_over_cycle(filename,legend):
    """Plot data from file vs line number in file. legend is list of strings to
    use for the plot legend, in the same order as the data in the file."""
    infile = open(filename)
    data = []
    for line in infile:
        values = line.split()
        data.append([float(val) for val in values])
    infile.close()
    fig,ax = plt.subplots()
    ax.plot(data,'.-')
    ax.set_xlabel("MCS")
    ax.legend(legend)
    return([fig,ax])

def plot_over_col1(filename,legend):
    """Plot data from file using the first column as the x-axis. legend is list
    of strings to use for the plot legend, in the same order as the data in the
    file."""
    infile = open(filename)
    infile.readline()
    x_data = []
    y_data = []
    for line in infile:
        values = line.split()
        x_data.append(float(values[0]))
        y_data.append([float(val) for val in values[1:]])
    infile.close()
    fig,ax = plt.subplots()
    ax.plot(x_data,y_data,'.')
    ax.set_xlabel("T")
    ax.legend(legend)
    return(np.array(y_data)[np.array(x_data).argsort()],[fig,ax])





#fig,ax = plt.subplots()
# Plot over T
#x_data = [Ti for is_good,Ti in zip([i==1e5 for i in MCS],T) if is_good]
#y_data = [Ei for is_good,Ei in zip([i==1e5 for i in MCS],expvals) if is_good]
# Plot over MCS
#x_data = [i for is_good,i in zip([Ti==1 for Ti in T],MCS) if is_good]
#y_data = [Ei for is_good,Ei in zip([Ti==1 for Ti in T],expvals) if is_good]
# Plot over MCS and T
#temps = np.sort(np.array(list(set(T))))

#for temp in temps:
#     x_data = [i for is_good,i in zip([Ti==temp for Ti in T],MCS) if is_good]
#     y_data = [Ei for is_good,Ei in zip([Ti==temp for Ti in T],expvals) if is_good]
#     y_data = [ai for is_good,ai in zip([Ti==temp for Ti in T],N_accepted_configs) if is_good]

     #x_data = np.sort(np.array(x_data)) # sort by MCS
     #y_data = [y for (x,y) in sorted(zip(x_data,y_data), key=lambda pair: pair[0])]
    # ax.semilogx(x_data,np.array(y_data),'.',label=str(temp))
#ax.legend()

def read_not_para_file():
    infile = open(outfilename)
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
    values = [E,M,N_accepted_configs]
    params = [L,MCS,T,spin0]
    return(values,params)

def energy_histogram(start):
    values,params = read_not_para_file()
    E = values[0][int(start):]
    MCS = params[1]-int(start)
    N = params[0]**2

    exp_E = np.sum(np.array(E))/MCS
    exp_E2 = np.sum(np.array(E)**2)/MCS

    Nbins = int((max(E)-min(E))/4) # bin width of 4J
    result = plt.hist(E,Nbins)
    plt.xlim((min(E),max(E)))

    sigmaE2 = exp_E2 - exp_E**2
    sigmaE = np.sqrt(sigmaE2)
    x = np.linspace(min(E),max(E),100)
    scale = len(E)*(result[1][1] - result[1][0])
    print(exp_E,sigmaE)
    plt.plot(x,stats.norm.pdf(x,exp_E,sigmaE)*scale)

def functions_of_MCS():
    [E,M,N_accepted_configs],params = read_not_para_file()
    MCS = np.linspace(1,params[1],int(params[1]))
    N = params[0]**2

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
        step = max([int(idxf/100),1])
        X += list(MCS[idx0:idxf:step])
        YE += list(exp_E[idx0:idxf:step]/N)
        YM += list(exp_M[idx0:idxf:step]/N)
        YN += list(np.array(N_accepted_configs)[idx0:idxf:step]/N)
        idx0 = idxf

    plt.semilogx(X,YE,'.-',label=r"$\langle E\rangle/N$")
    plt.semilogx(X,YM,'.-',label=r"$\langle |M|\rangle/N$")
    plt.semilogx(X,YN,'.-',label="accepted configs")
    plt.xlabel("Monte Carlo cycles")
    # plot every 1000th value from 1000 to
    plt.legend()




#energy_histogram(100000)
#functions_of_MCS()
#ax.plot(x_data,y_data,'.')
#ax.legend(["E","M","$E^2$","$M^2$","$C_V$","$\chi$"])


#[fig,ax] = plot_over_cycle("10x10_mcs.txt",["E","M"])
y_data,[fig,ax] = plot_over_col1(outfilename,["E","M","$E^2$","$M^2$","$C_V$","$\chi$"])

T = np.linspace(0.1,3,30)
exact = Ising_2x2.Ising_2x2(T)

difference = -np.array([exact.E/4,exact.abs_M/4,exact.C_V/4,exact.chi/4]).transpose() + np.delete(np.delete(y_data,2,1),2,1)

plt.gca().set_prop_cycle(None)
exact.plot_2x2(ax)
ax.legend()

plt.figure()
plt.plot(T,difference,'-')#np.array([exact.E/4,exact.abs_M/4,exact.C_V/4,exact.chi/4]).transpose(),'.-')
plt.legend([r'$\langle E\rangle/N$',r'$\langle |M|\rangle/N$',r'$\langle C_V\rangle/N$',r'$\langle \chi\rangle/N$'])
plt.axhline(0,0,1,ls='--',c='k')
plt.xlabel('T')
plt.ylabel('deviation from theory')

plt.show()
