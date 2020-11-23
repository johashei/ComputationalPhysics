all: compile link

compile:
	g++-10 -fopenmp -O3 -c main.cpp Ising.cpp not_para.cpp -std=c++11

link:
	g++-10 -fopenmp -O3 -o main_omp main.o Ising.o
	g++-10 -O3 -o main_notpara not_para.o Ising.o
