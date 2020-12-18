The files ``SIRS.cpp`` and ``SIRS.hpp`` define the class SIRS, which perform Monte Carlo simulation of the SIRS model. Descriptions of the methods, as well as how to expand the class are given in the program file. The class needs to be compiled with  `C++11`. 

---

The class SIRS defined in ``SIRS.py`` simulates the SIRS model with Runge-Kutta 4, as implemented in ``ODESolver.py``. More information about this class can be found in ([Langtanged and Sundnes, 2017](https://www.uio.no/studier/emner/matnat/ifi/IN1900/h17/ressurser/slides/ode\_systems\_short.pdf)).

---

The program files ``main.cpp`` and ``main.py`` were used when creating the figures for the report. Running `main.py` also compiles, links and runs `main.cpp` with the same arguments. Note that `main.cpp` is linked to `SIRS.o`, so the class code `SIRS.cpp` should be compiled beforehand.  

They take 13 command line arguments:

   - `a0` : float, rate of transmission $\beta$
   - `b` : float, rate of recovery $\gamma$  
   - `c` : foat, rate of immunity loss $\lambda$
   - `S`,`I`,`R` : int, initial values for $S,I,R$
   - `Tfinal` : float, simulation time
   - `d` : float, rate of death $\mu$ unrelated to the disease
   - `dI` : float, rate of death $\mu_I$ due to the disease
   - `e` : float, rate of birth $\nu$
   - `A` : float, Amplitude of the seasonal variation
   - `omega` : float, angular frequency of the seasonal variation
   - `f0` : float, a variable used in the vaccination rate function. 

Both programs define a function `f`, which describes the vaccination rate as a function of time. When modifying this function, care should be taken that it is equivalent in both programs.     

The SIRS class called by `main.cpp` takes an optional seed value for its random number generator in its initializer. This value needs to be specified in the `main.cpp` code.  Using the variable `seed` gives a different seed each time the code in run.

The program `main.py` reads the file created by the SIRS class in `main.cpp` and plots the ODE and Monte Carlo solutions together in the same plot.  

---








