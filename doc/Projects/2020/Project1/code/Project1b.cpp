#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "time.h"

using namespace std;

// Declarations of variables
clock_t start, finish; // declare start and final time

// Declarations of functions
double f(double); //  right hand side function
double v(double); //  analytical solution function
double error(int,double*,double*); // maximum error
double* general_algo(int,double); // general tridiagonal solver
double* specific_algo(int,double); // 2nd order DE solver


int main(int argc, char *argv[]){

  // read filename and number of steps from command line
  string outfilename;
  int n;
  if( argc <= 3){
    cout << "Bad usage. Use " << argv[0] << " <output file name> <number of steps> <algorithm>" << endl;
    exit(1);
  }
  outfilename = argv[1];
  n = atoi(argv[2]);
  double h = 1.0/(n+1); // stepsize x_0 = 0 and x_n+1 = 1

  // Run and time algorithm:
  double* u;
  if(strcmp(argv[3], "general") == 0){
    start = clock(); // start timer
    u = general_algo(n,h);
  }
  else if(strcmp(argv[3], "specific") == 0){
    start = clock(); // start timer
    u = specific_algo(n,h);
  }
  else{
    cout << "Bad usage. <algorithm> should be 'general' or 'specific'." << endl;
    exit(1);
  }
  finish = clock(); // stop timer
  double Time = (double(finish - start)/CLOCKS_PER_SEC);

  // Analytic solution
  double *v_array = new double[n+2];
  for(int i=0 ; i<n+2 ; i++){
    v_array[i] = v(i*h);
  }

  // Maximum error
  double maxerror = error(n,u,v_array);

  // Write results to file
  ofstream ofile;
  ofile.open(outfilename);
  ofile << "Solution to differential equation with "<< argv[3] << " algorithm:" << endl;
  ofile << n << endl;
  for(int i=0; i<n+2; i++){
    ofile << setw(15) << setprecision(8) << u[i] << " ";
    ofile << setw(15) << setprecision(8) << v_array[i] << endl;
  }
  ofile << "Max relative error and log10(h):" << endl;
  ofile << setw(15) << setprecision(8) << maxerror << " ";
  ofile << setw(15) << setprecision(8) << log10(h) << endl;
  ofile << "Time taken by the algorithm [s]:" << endl;
  ofile << Time << endl;
  ofile.close();

  delete [] u; delete [] v_array;

  return 0;
}

double* general_algo(int n, double h){
  // Creating arrays:
  double *a = new double [n-1];
  double *b = new double [n];
  double *c = new double [n-1];
  double *g = new double [n];
  double *u = new double [n+2];
  double *btilde = new double [n];
  double *gtilde = new double [n];

  // Initializing:
  for (int i = 0; i<n-1; i++){
    a[i] = -1;
    b[i] = 2;
    c[i] = -1;
  }
  b[n-1] = 2;

  for(int i=0 ; i<n+1; i++){
    g[i] = h*h*f(i*h);
  }
  // Solution, assuming boundary u(0) = u(1) = 0
  u[0] = 0;
  u[n+1] = 0;

  btilde[0] = b[0];
  gtilde[0] = g[0];

  // Forward substitution -- 5N FLOPs
  for(int i=1 ; i<n ; i++){
    double row_reduction_factor = a[i-1]/btilde[i-1]; // used twice so saves one FLOP
    btilde[i] = b[i] - c[i-1]* row_reduction_factor;
    gtilde[i] = g[i] - gtilde[i-1]* row_reduction_factor;
    cout << btilde[i] << " ";
  }
  cout << endl;
  //Backward substitution -- 3N FLOPs
  for(int i=n+1 ; i>1 ; i--){
    u[i-1] = (gtilde[i-2] - c[i-2]*u[i])/btilde[i-2]; // -2 because they start at 0 not 1
  }
  return u;
}

double* specific_algo(int n, double h){
  // Something wrong somewhere see file for n=10 ?
  // Special case where a_i = c_i = -1 and b_i = 2 for all i
  // Creating arrays:
  double *g = new double [n];
  double *u = new double [n+2];
  double *edlitb = new double [n]; // 1/btilde, because * is faster than /
  double *gtilde = new double [n];

  // Initializing
  for(int i=0 ; i<n+1; i++){
    g[i] = h*h*f(i*h);
  }
  // Solution, assuming boundary u(0) = u(1) = 0
  u[0] = 0;
  u[n+1] = 0;

  gtilde[0] = g[0];

  // Precalculating 1/btilde
  edlitb[0] = 1.0/2;
  for(int i=1 ; i<n ; i++){
    edlitb[i] = double(i+1)/(i+2);
    cout << 1/edlitb[i] << " ";
  }
  cout << endl;
  // Forward substitution -- 2N FLOPs
  for(int i=1 ; i<n ; i++){
    gtilde[i] = g[i] + gtilde[i-1]*edlitb[i-1];
  }
  // Backward substitution -- 2N FLOPs
  for(int i=n+1 ; i>1 ; i--){
    u[i-1] = (gtilde[i-2] + u[i])*edlitb[i-2];
  }
  delete [] g; delete [] gtilde; delete [] edlitb;
  return u;
}

double error(int n, double *u, double *v){
  // find the maximum relative error
  double maxerror = -1e4; // probably smaller than any max error
  for(int i=1 ; i<n+1 ; i++){
    double epsilon = log10(abs((u[i] - v[i])/v[i])); // Why log10???
    //cout << epsilon << " ";
    if(epsilon >= maxerror){
      maxerror = epsilon;
    }
  }
  return maxerror;
}

double f(double x){
  return 100*exp(-10*x);
}

double v(double x){
  return 1 - (1-exp(-10))*x - exp(-10*x);
}
