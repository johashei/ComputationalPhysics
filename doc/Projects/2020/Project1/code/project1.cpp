#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include "time.h"
#include "lib.h"

using namespace std;

// Declarations of variables
clock_t start, finish; // declare start and final time

// Declarations of functions
double f(double); //  right hand side function
double v(double); //  analytical solution function
double error(int,double*,double*); // maximum error
double* general_algo(int,double); // general tridiagonal solver
double* specific_algo(int,double); // 2nd order DE solver
double* LU_algo(int,double); // LU decomp and solver from lib.h


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
    cout << "Solving with general algo ..." << flush;
    start = clock(); // start timer
    u = general_algo(n,h);
    cout << " Done." << endl;
  }
  else if(strcmp(argv[3], "specific") == 0){
    cout << "Solving with specific algo ..." << flush;
    start = clock(); // start timer
    u = specific_algo(n,h);
    cout << " Done." << endl;
  }
  else if(strcmp(argv[3], "LU") == 0){
    cout << "Solving with LU decomposition algo ..." << flush;
    start = clock(); // start timer
    u = LU_algo(n,h);
    cout << " Done." << endl;
  }
  else{
    cout << "Bad usage. <algorithm> should be 'general', 'specific' or 'LU'." << endl;
    exit(1);
  }
  finish = clock(); // stop timer
  //cout << finish << " "<< start << " " << CLOCKS_PER_SEC << endl;
  double Time = (double(finish - start)/CLOCKS_PER_SEC)*1e3; // time in ms


  // Analytic solution
  cout << "Calculating analytic solution ..." << flush;
  double *v_array = new double[n+2];
  for(int i=0 ; i<n+2 ; i++){
    v_array[i] = v(i*h);
  }
  cout << " Done." << endl;

  // Maximum error
  cout << "Calculating maximum relative error ..." << flush;
  double maxerror = error(n,u,v_array);
  cout << " Done." << endl;

  if(outfilename != "None"){// Write results to file
    //start = clock(); // New clock to time write to file
    cout << "Writing results to file ..." << flush;
    ofstream ofile;
    ofile.open(outfilename);
    ofile << "Solution to differential equation with "<< argv[3] << " algorithm:\n";
    ofile << n << "\n";
    for(int i=0; i<n+2; i++){
    ofile << setw(15) << setprecision(8) << u[i] << " ";
    ofile << setw(15) << setprecision(8) << v_array[i] << "\n";
  }
    ofile << "Max relative error and log10(h):\n";
    ofile << setw(15) << setprecision(8) << maxerror << " ";
    ofile << setw(15) << setprecision(8) << log10(h) << "\n";
    ofile << "Time taken by the algorithm [ms]:\n";
    ofile << Time << "\n";
    ofile.close();
    cout << " Done." << endl;
    //finish=clock();
    //cout << double(finish - start)/CLOCKS_PER_SEC << endl;
  }
  else{// Display Time. Used for timing experiment { } symbols used to mark position of time value
    cout << "Time taken by the algorithm [ms]: {" << Time << "}" << endl;
  }

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
  // Initializing:
  for(int i=0 ; i<n+1; i++){
    g[i] = h*h*f(i*h);
  }
  for (int i = 0; i<n-1; i++){
    a[i] = -1;
    b[i] = 2;
    c[i] = -1;
  }
  b[n-1] = 2;
  // Solution, assuming boundary u(0) = u(1) = 0
  u[0] = 0;
  u[n+1] = 0;
  // Forward substitution -- 5n FLOPs; 8n fetches; 3n writes
  for(int i=1 ; i<n ; i++){
    double row_reduction_factor = a[i-1]/b[i-1]; // used twice so saves one FLOP
    b[i] = b[i] - c[i-1]* row_reduction_factor;
    g[i] = g[i] - g[i-1]* row_reduction_factor;
  }
  //Backward substitution -- 3n FLOPs; 4n fetches; 1n writes
  for(int i=n+1 ; i>1 ; i--){
    u[i-1] = (g[i-2] - c[i-2]*u[i])/b[i-2]; // -2 because they start at 0 not 1
  }
  delete [] a; delete [] b; delete [] c; delete [] g;
  return u;
}

double* specific_algo(int n, double h){
  // Special case where a_i = c_i = -1 and b_i = 2 for all i
  // Creating arrays:
  double *g = new double [n];
  double *u = new double [n+2];
  double *edlitb = new double [n]; // 1/btilde, because * is faster than /
  // Initializing
  for(int i=0 ; i<n+1; i++){
    g[i] = h*h*f(i*h);
  }
  // Solution, assuming boundary u(0) = u(1) = 0
  u[0] = 0;
  u[n+1] = 0;
  // Precalculating 1/btilde
  edlitb[0] = 1.0/2;
  for(int i=1 ; i<n ; i++){
    edlitb[i] = double(i+1)/(i+2);
  }
  // Forward substitution -- 2n FLOPs; 3n fetches; 1n writes
  for(int i=1 ; i<n ; i++){
    g[i] = g[i] + g[i-1]*edlitb[i-1];
  }
  // Backward substitution -- 2n FLOPs; 3n fetches; 1n writes
  for(int i=n+1 ; i>1 ; i--){
    u[i-1] = (g[i-2] + u[i])*edlitb[i-2];
  }
  delete [] g; delete [] edlitb;
  return u;
}

double* LU_algo(int n, double h){
  // Creating arrays:
  double **A;
  double *g = new double[n];
  double *u = new double[n+2];
  // Allocate space for lib function args
  int *indx = new int[n];
  double d;
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
  // Initialise right hand side g:
  for(int i=0 ; i<n; i++){
    g[i] = h*h*f(i*h);
  }
  // Call functions from lib
  ludcmp(A,n,indx,&d); // LU decomposition
  lubksb(A,n,indx,g); // Solves Au = g, stores u in g
  // Put solution in array u
  u[0] = 0;
  u[n+2] = 0;
  for(int i=0; i<n; i++){
    u[i+1] = g[i];
  }
  // Delete matrix A
  for(int i=0; i<n; i++){
    delete [] A[i] ;
  }
  delete [] A;
  delete [] g;

  return u;
}

double error(int n, double *u, double *v){
  // find the maximum relative error
  double maxerror = log10(abs((u[1] - v[1])/v[1])); // first error
  for(int i=2 ; i<n+1 ; i++){
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
