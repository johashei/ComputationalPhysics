#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "time.h"
#include "lib.h"

clock_t start, finish; // declare start and final times
double **A; // declare tridiagonal matrix
double *g; // declare pointer to right hand side then solution

double f(double); // declare right hand side function

int main(int argc, char *argv[]){
  // read filename and number of steps from command line
  string outfilename;
  int n;
  if( argc <= 2){
    cout << "Bad usage. Use " << argv[0] << " <output file name> <number of steps>" << endl;
    exit(1);
  }
  outfilename = argv[1];
  n = atoi(argv[2]);
  double h = 1.0/(n+1); // stepsize x_0 = 0 and x_n+1 = 1

  if(n>=1e4){
    cout << "order of n should not exceed 10^3." << endl;
    exit(1);
  }

  // Initialise matrix A :
  A = new double*[n];
  for(int i=0; i<n; i++){
    A[i] = new double[n];
    for(int j=0; j<n; j++){
      A[i][j] = 0.0;
    }
  }
  // Fill in tridiagonal elements :
  for(int i=0; i<n-1; i++){
    A[i][i] = 2;
    A[i][i+1] = -1;
    A[i+1][i] = -1;
  }
  A[n-1][n-1] = 2;

  // Initialise right hand side b:
  g = new double[n];
  for(int i=0 ; i<n+1; i++){
    g[i] = h*h*f(i*h);
  }

  // Allocate space for lib function args
  int *indx = new int[n];
  double d;

  start = clock();
  // Call functions from lib
  ludcmp(A,n,indx,&d); // LU decomposition
  lubksb(A,n,indx,g); // Solves Au = g, stores u in g
  finish = clock();
  double Time = (double(finish - start)/CLOCKS_PER_SEC);


  // Write results to file
  ofstream ofile;
  ofile.open(outfilename);
  ofile << "Solution to differential equation with LU decomposition from lib.cpp \n:";
  ofile << n << "\n";
  ofile << 0 << "\n";
  for(int i=0; i<n; i++){
    ofile << setw(15) << setprecision(8) << g[i] << "\n";
  }
  ofile << 0 << "\n";
  ofile << "Time taken by the algorithm [s]:\n";
  ofile << Time << "\n";
  ofile.close();

  delete [] g;
  for(int i=0; i<n; i++){
    delete [] A[i] ;
  }
  delete [] A;

  return 0;
}

double f(double x){
  return 100*exp(-10*x);
}
