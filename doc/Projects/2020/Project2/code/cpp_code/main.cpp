#include "EigenValueSolver.hpp"
#include <cmath>
#include <string>
#include <ctime>
#include <armadillo>

using namespace std;
using namespace arma;

dmat tridiagonal_Toeplitz_matrix_sym(int N, double offd, double d);

double potential(double rho);

int main(int argc, char const* argv[]){

  clock_t start, end;
  double algotime;
  uword eigval_idx_max = 3; // Display the four lowest eigenvalues

  // Read args from command line:
  double rho_0 = atoi(argv[1]);
  double rho_N = atoi(argv[2]); // Approximation of rho_max = inf
  int N = atoi(argv[3]);
  string algorithm = string(argv[4]);
  uword M; // Must be declared outside if test to avoid compiler error
  if(algorithm=="Lanczos"){
    M = atoi(argv[5]); 
  }

  double h = (rho_N-rho_0)/N;

  mat V = zeros<mat>(N-1,N-1); // Diagonal matrix with potentials for each rho
  for(int i=0; i<N-1; i++){
    V(i,i) = potential(rho_0 + (i+1)*h);// Fill potential matrix
  }

  double a = -1.0/(h*h);
  double d = 2.0/(h*h);
  dmat A = tridiagonal_Toeplitz_matrix_sym(N,a,d) + V;


  // Solve with Armadillo, Lapack:{
  cx_vec eigval;
  cx_mat eigvec;
  start = clock();
  eig_gen(eigval,eigvec,A);
  end = clock();
  algotime = (double) (end-start)/CLOCKS_PER_SEC;
  vec Reigval = sort(real(eigval)); // take real part and sort
  cout << "\n With Lapack:"<<endl;
  cout << real(Reigval.subvec(0,eigval_idx_max))<<endl;
  cout << "Time used: " << algotime << " s\n" << endl;
  //cout << real(eigvec) << endl;
  //}

  if(algorithm == "Jacobi"){
    cout << "Jacobi's algorithm:" << endl;
    double tolerance = 1e-8;
    int maxiter = 1e6;
    Jacobi_rotation jacobi_solver;
    jacobi_solver.init(A,tolerance,maxiter);
    start = clock();
    jacobi_solver.solve();
    end = clock();
    algotime = (double) (end-start)/CLOCKS_PER_SEC;

    cout << "Eigenvalues:" << endl;
    cout << jacobi_solver.Get_eigvals(0,eigval_idx_max) << endl;
    cout << "Time used: " << algotime << " s" << endl;
  }

  if(algorithm == "Lanczos"){
    cout << "Lanczos' algorithm:" << endl;
    Lanczos lanczos_solver;
    lanczos_solver.init(A,M);
    start = clock();
    lanczos_solver.solve();
    end = clock();
    algotime = (double) (end-start)/CLOCKS_PER_SEC;

    cout << "Eigenvalues:" << endl;
    cout << lanczos_solver.Get_eigvals(0,eigval_idx_max) << endl;
    cout << "Time used: " << algotime << " s" << endl;
  }

  return 0;
}

dmat tridiagonal_Toeplitz_matrix_sym(int N, double offd, double d){
  vec A1 = zeros<vec>(N-1);
  A1(0) = d;
  A1(1) = offd;
  dmat A = toeplitz(A1);
  return A;
}

double potential(double rho){
  double V = rho*rho;
  return V;
}
