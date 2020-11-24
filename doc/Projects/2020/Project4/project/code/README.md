The code can be compiled and linked using ```make all```.

---

```main.cpp``` is a parallelized program for solving the two dimentional Ising model. The number of threads ```MAX_THREADS``` is set whithin the program file. The parameters are supplied as command line arguments, in the following order:

```int L```: side length of lattice

 ```int MCS``` : number of Monte Carlo cycles

```int equilibration``` : number of Monte Carlo cycles to discard when computing expectation values

```double T_start``` : first temperature

```double T_end``` : last temperature

```double T_step``` : temperature step

```string spin0``` : initial state of the lattice should be "aligned", "random", or "default". Default is aligned if T<1.5 and random otherwise

```string outfilename``` : name of file to generate 

The first line of the generated file displays the input parameters

 L	MCS	equilibration	spin0

The remaining lines display the various expectation values per spin for each temperature

T 	E 	|M|	 E<sup>2</sup>   M<sup>2</sup> 	C<sub>V</sub> 	𝜒

---

The program ```not_para.cpp``` is not parallelized, and simulates the Ising model for a single temperature T. The command line arguments are ```L```,```MCS```,```T```,```spin0``` and ```outfilename```. 

The first line in the output file displays the four first arguments. The remaining lines display the energy, magnetization and number of accepted configurations for each Monte Carlo cycle.

---

Both of these programs use the ```Ising_2d``` class. See the source file ```Ising.cpp``` for information on this class. 

---

The python class ```Ising_2x2``` calculates the analytical expectation values for the 2x2 spin lattice. 

---

The calls used to generate the data in the report are listed in ```commands.txt```. This file also contains the terminal output of the runs, including the runtime for each call. 

The figures  were created using the programs ```mcplots.py```,```histograms.py```, ```mainplots.py``` and ```interface.py```. These were made specifically for the data in this report and may need to be modified to work with other datafiles. 

---

All programs were compiled and run using gcc version 10.2.0 (Homebrew GCC 10.2.0) on a 4-core Intel Core i7 MacBook Pro with macOS 10.15.4

