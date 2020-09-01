#include "general_tridiagonal.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

void general_tridiagonal::Initialize(int N, double * a, double * b, double * c, double f(double x)){
  m_N = N;
  // Setting up arrays
  m_a = new double [N-1]; m_b = new double [N]; m_c = new double [N-1];
  m_g = new double [N]; m_u = new double [N+2];
  // Matrix elements
  for(int i = 0; i<m_N-1; i++){
    m_a[i] = a[i];
    m_b[i] = b[i];
    m_c[i] = c[i];
  }
  m_b[m_N-1] = b[m_N-1];
  // Right hand side vector, assuming function limits 0, 1
  double h = 1.0/(m_N+1); // stepsize x_0 = 0 and x_n+1 = 1
  for(int i=0 ; i<m_N; i++){
    m_g[i] = h*h*f(i*h);
  }
  // Solution, assuming boundary u(0) = u(1) = 0
  m_u[0] = 0;
  m_u[m_N+1] = 0;
}

void general_tridiagonal::Solve(){

  double * btilde = new double [m_N]; // function variable
  double * gtilde = new double [m_N]; // function variable
  btilde[0] = m_b[0];
  gtilde[0] = m_g[0];

  // Forward substitution -- 5N FLOPs
  for(int i=1 ; i<m_N ; i++){
    double row_reduction_factor = m_a[i-1]/btilde[i-1]; // used twice so saves one FLOP
    btilde[i] = m_b[i] - m_c[i-1]* row_reduction_factor;
    gtilde[i] = m_g[i] - gtilde[i-1]* row_reduction_factor;
  }
  //Backward substitution -- 3N FLOPs
  for(int i=m_N+1 ; i>1 ; i--){
    m_u[i-1] = (gtilde[i-1] - m_c[i-2]*m_u[i])/btilde[i-2];
  }
  delete [] btilde; delete [] gtilde;
}

void general_tridiagonal::Writefile(string outfilename){
  ofstream ofile;
  ofile.open("./1b_data/"+outfilename);
  ofile << "Solution to differential equation with general tridiagonal algorithm:" << endl;
  ofile << m_N << endl;
  for(int i=0; i<m_N+2; i++){
    ofile << setw(15) << setprecision(8) << m_u[i] << endl;
  }
}
