import matplotlib.pyplot as plt
import numpy as np
import sys,os,subprocess
import Ising_2x2

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
    #ax.legend(legend)
    return([fig,ax])



#[fig,ax] = plot_over_cycle("2x2EM.txt",["E","M"])
[fig,ax] = plot_over_col1("2x2_exp.txt",["E","M","$E^2$","$M^2$","$C_V$","$\chi$"])

T = np.linspace(0.1,3)
example = Ising_2x2.Ising_2x2(T)
plt.gca().set_prop_cycle(None)
example.plot_2x2(ax)

ax.legend()
plt.show()
