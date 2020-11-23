main.cpp is a parallelized program for solving the two dimentional Ising model. The number of threads MAX_THREADS is set whithin the program file. The parameters are supplied as command line arguments, in the following order:

```
int L : side length of lattice
int MCS : number of Monte Carlo cycles
int equilibration : number of Monte Carlo cycles to discard when computing expectation values
double T_start : first temperature
double T_end : last temperature
double T_step : temperature step
string spin0 : initial state of the lattice should be "aligned", "random", or "default". Default is aligned if T<1.5 and random otherwise
string outfilename : name of file to generate 
```

