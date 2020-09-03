all: compile execute

compile:
	c++ -c -Wall project1.cpp general_tridiagonal.cpp
	c++ -o Project1 project1.o general_tridiagonal.o
	c++ -o Project1b Project1b.cpp

execute:
	./Project1


all_alternative: compile_alternative execute

compile_alternative:
	c++ -o Project $(wildcard *.cpp)

reset:
	rm ./P1data/*
