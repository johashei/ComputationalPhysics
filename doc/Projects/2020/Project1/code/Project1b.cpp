#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

double f(double); // declare right hand side function
double v(double); // declare analytical solution function
double error(int,double*,double*); // declare function to calculate maximum error

// Declare the two different algorithms
double* general_algo(int,double);
double* specific_algo(int,double,double,double,double);

int main(int argc, char *argv[]){

  // read filename and number of steps from command line
  string outfilename;
  int n;
  if( argc <= 2){
    cout << "Bad usage. Use " << argv[0] << " <output file name>  <number of steps>" << endl;
    exit(1);
  }
  outfilename = argv[1];
  n = atoi(argv[2]);
  double h = 1.0/(n+1); // stepsize x_0 = 0 and x_n+1 = 1

  // Algorithm:
  double* u;
  u = general_algo(n,h);
  // u = specific_algo(n,h,-1,2,-1);

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
  ofile << "Solution to differential equation with general tridiagonal algorithm:" << endl;
  ofile << n << endl;
  for(int i=0; i<n+2; i++){
    ofile << setw(15) << setprecision(8) << u[i] << " ";
    ofile << setw(15) << setprecision(8) << v_array[i] << endl;
  }
  ofile << "Max relative error and log10(h):" << endl;
  ofile << setw(15) << setprecision(8) << maxerror << " ";
  ofile << setw(15) << setprecision(8) << log10(h) << endl;
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
  b[n-1] = -2;

  for(int i=0 ; i<n+1; i++){
    //cout << i*h << endl;
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
  }
  //Backward substitution -- 3N FLOPs
  for(int i=n+1 ; i>1 ; i--){
    u[i-1] = (gtilde[i-2] - c[i-2]*u[i])/btilde[i-2];
  }

  //delete[] a; delete[] b; delete[] c; delete[] g;
  //delete[] btilde; delete[] gtilde;

  return u;
}

double* specific_algo(int n, double h, double a, double b, double c){
  // Creating arrays:
  double *g = new double [n];
  double *u = new double [n+2];
  double *btilde = new double [n];
  double *gtilde = new double [n];

  // Initializing

  for(int i=0 ; i<n+1; i++){
    //cout << i*h << endl;
    g[i] = h*h*f(i*h);
  }
  // Solution, assuming boundary u(0) = u(1) = 0
  u[0] = 0;
  u[n+1] = 0;

  btilde[0] = b;
  gtilde[0] = g[0];

  // Precalculating elements:
  // btilte_i = b - c*a/btilde_i-1 =
  return 0;

}

double error(int n, double *u, double *v){
  // find the maximum relative error
  double maxerror = 0;
  for(int i=0 ; i<n+2 ; i++){
    double epsilon = log10(abs((u[i] - v[i])/v[i]));
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
