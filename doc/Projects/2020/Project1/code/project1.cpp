#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "general_tridiagonal.hpp"

using namespace std;

double f(double x);

int main(int argc, char *argv[]) {
  // read filename and number of steps from command line
  string outfilename;
  int N;
  if( argc <= 2){
    cout << "Bad usage. Use " << argv[0] << " <output file name>  <number of steps>" << endl;
    exit(1);
  }
  outfilename = argv[1];
  N = atoi(argv[2]);

  // General case
  double * a = new double [N-1];
  double * b = new double [N];
  double * c = new double [N-1];

  for (int i = 0; i<N-1; i++){
    a[i] = -1;
    b[i] = 2;
    c[i] = -1;
  }
  b[N] = -2;


  general_tridiagonal solver; // Create solver of type general_tridiagonal
  solver.Initialize(N,a,b,c,f);
  solver.Solve();
  solver.Writefile(outfilename);

  delete [] a; delete [] b; delete [] c;
  //cout << "done" << endl;
  return 0;
}

double f(double x){
  return 100*exp(-10*x);
}
