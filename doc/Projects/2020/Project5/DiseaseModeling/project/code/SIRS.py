import matplotlib.pyplot as plt
import numpy as np
import ODESolver

if __name__!="__main__":
    print("SIRS.py imported")

class SIRS(object):
    """Solve the SIRS model as a system of coupled ODEs using the RK4 algorithm"""

    def __init__(self,S,I,R):
        """Class initializer: constructs the state of the system at time t=0.
        Arguments:  S : Initial number of susceptible
                    I : Initial number of infected
                    R : Initial number of recovered
        """
        self.state = [S,I,R]

    def equations(self,state,t):
        """Create the system of ODEs to be used by ODESolver
        """
        S,I,R = state
        N = S+I+R
        for transition in self.transitions:
            exec(transition)

        DS = eval(self.DS_string)
        DI = eval(self.DI_string)
        DR = eval(self.DR_string)
        return [DS,DI,DR]

    def Solve_ODE(self,final_time,dt,terminate=None):
        """Use ODESolver's RungeKutta4 method to solve the system of equations
        Arguments:  final_time : length of the simulation
                    dt : time step
                    terminate (optional) : condition upon which to break the solver loop
        """
        n = round(final_time/dt + 1) # n must be an integer
        print(dt,n)
        t = np.linspace(0,final_time,n)
        ODE = ODESolver.RungeKutta4(self.equations)
        ODE.set_initial_condition(self.state)
        SIR,t = ODE.solve(t,terminate)

        return(SIR,t)

    def Basic_SIRS(self, rStoI,rItoR,rRtoS):
        """Define the equations of the basic SIRS model. This shoul be called
        first after the constructor.
        Arguments:  rStoI : transmission parameter
                    rItoR : recovery parameter
                    rRtoS : immunity loss parameter
        """
        self.a0 = rStoI
        self.b = rItoR
        self.c = rRtoS
        self.time = 0
        self.a = lambda t: self.a0 # Without vital dynamics
        # ROOT-like and ugly. Find a better way to do this.
        self.transitions = ["StoI = self.a(t)*S*I/N",
                            "ItoR = self.b*I",
                            "RtoS = self.c*R"]
        # These will be evaluated in the method "equations"
        self.DS_string = "RtoS - StoI"
        self.DI_string = "StoI - ItoR"
        self.DR_string = "ItoR - RtoS"

    def Vital_dynamics(self,rBirth,rDeath,rIDeath):
        """Add the terms for vital dynamics to the equations
        Arguments:  rBirth : birth rate
                    rDeath : rate of death not related to the disease
                    rIDeath : rate of death caused by the disease
        """
        self.e = rBirth
        self.d = rDeath
        self.dI = rIDeath
        self.DS_string += "+ self.e*(S+I+R) - self.d*S"
        self.DI_string += "- self.d*I - self.dI*I"
        self.DR_string += "- self.d*R"

    def Seasonal_variation(self,amplitude, omega):
        """Change the parameter self.a to account for seasonal variation
        Arguments:  amplitude : maximum deviation from the mean transmission rate
                    omega : angular frequency of the oscillation
        """
        self.A = amplitude
        self.w = omega
        self.a = lambda t: self.A * np.cos(self.w*t) + self.a0

    def Vaccination(self,rate):
        """Add the term for vaccination to the equations
        Argument:   rate(t)  : function of time describing the vaccination rate.
        """
        self.rVacc = lambda t: rate(t)
        self.transitions.append("StoR = self.rVacc(t)*(S>0)")
        self.DS_string += "- StoR"
        self.DR_string += "+ StoR"
