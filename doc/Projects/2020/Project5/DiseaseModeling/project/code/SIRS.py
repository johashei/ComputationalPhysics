"""Not part of the project. Made to test results of project 5 in FYS4150, 2020"""

import matplotlib.pyplot as plt
import numpy as np
import ODESolver

print("SIRS.py imported")

def SIRSSolver(S_0,I_0,R_0, beta,b,c, dt,t_i,t_f, terminate=None):
    '''On the time interval [t_i,t_f]: Solve the differential equations DS(x,t),
    DI(x,t) and DR(x,t) in the SIRS model. Return S(t), I(t), R(t), t.'''

    n = round((t_f - t_i)/dt + 1) # n must be an integer
    t = np.linspace(t_i, t_f, n)
    u0 = [S_0, I_0, R_0]
    def f(u,t):
        S, I, R = u
        DS = - beta * S * I + c*R
        DI = beta * S * I - b * I
        DR = b * I - c*R
        return [DS, DI, DR]

    # Using ODESolver:
    ODE = ODESolver.RungeKutta4(f)
    ODE.set_initial_condition(u0)
    SIR, t = ODE.solve(t, terminate)
    # Returns an array with [S(t), I(t), R(t)]arrays and an array with t.

    S, I, R = SIR[:,0], SIR[:,1], SIR[:,2]

    return S, I, R, t


def visualizer(S,I,R, t, pop_unit=None, t_unit=None):
    '''visualize S(t), I(t) and R(t) in one plot.'''
    unit = [pop_unit, t_unit]
    for i in [0,1]:
        if unit[i] is None:
            unit[i] = ''
        else:
            unit[i] = '(in {unit})'.format(unit=unit[i])
    pop_unit, t_unit = unit[0], unit[1]

    plt.plot(t,S,'xkcd:blue', t,I,'xkcd:vomit', t,R,'xkcd:hospital green')
    plt.legend(['Susceptibles', 'Infected', 'Recovered'])
    plt.xlabel('time %s' %t_unit)
    plt.ylabel('people %s' %pop_unit)
    # Not including plt.show() so the plots can be used in the main program.


if __name__ == '__main__':
    # Setting Variables:
    S_0 = 300 ; I_0 = 100 ; R_0 = 0
    dt = 0.01 ; t_i = 0 ; t_f = 30
    b = 1 ; beta = 0.01 ; c = 0.5

    # Making plots:
    plt.figure('Simulation of the spreading of a disease by the SIRS model')

    plt.subplot()
    S, I, R, t = SIRSSolver(S_0,I_0,R_0, beta,b,c, dt,t_i,t_f)
    visualizer(S,I,R, t, t_unit = 'days')
    plt.title("$\\beta=${beta}".format(beta=beta))

    plt.subplots_adjust(left=0.12, bottom=0.10, right=0.99, top=0.95, hspace=0.40)
    plt.show()
