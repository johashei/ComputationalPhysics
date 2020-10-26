The python scrips <code>interface.py</code> compiles, links and runs the main program, reads the generated data file and plots the result. It requires all .o files apart from the main program to already exist. Running

<code> 
Terminal > python interface.py main </code></p>
<code>
Terminal > Input name of datafile : example.dat
</code>

will compile and run the program `main.cpp`, with "example.dat" as an argument, then read the file `example.dat` and create a plot `example.eps` from the data.

The programs `SolarSystem`, `TwobodySystem` and `simple` can be used as `main`. Note that the file name is not taken as an argument by `simple`, so use `simple.dat` or change it in the code.


The results were produced with Apple clang version 11.0.3 (clang-1103.0.32.59) in macOS 10.15.4
