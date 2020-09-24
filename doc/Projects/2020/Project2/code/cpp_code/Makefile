all: main test

main:
	c++ -o main.out -Wall EigenValueSolver.cpp Jacobi_rotation.cpp main.cpp -larmadillo

test:
	c++ -o test.out -Wall EigenValueSolver.cpp Jacobi_rotation.cpp test-functions.cpp tests-main.cpp -larmadillo -std=c++11
