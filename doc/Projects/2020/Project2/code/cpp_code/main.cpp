#include "EigenValueSolver.hpp"
#include <cmath>
#include <string>
#include <ctime>
#include <armadillo>

using namespace std;
using namespace arma;

dmat tridiagonal_Toeplitz_matrix(int N, double lower, double center, double upper);

int main(int argc, char const* argv[]){

  int N = atoi(argv[1]);
  string filename = string(argv[2]);

  double h = 1.0/N;
  double a = -1.0/(h*h);
  double d = 2.0/(h*h);
  dmat A = tridiagonal_Toeplitz_matrix(N,a,d,a);
  cout << A << endl;

  double tolerance = 1e-8;
  int maxiter = 1e3;
  Jacobi_rotation buckling_beam_solver;
  buckling_beam_solver.init(A,tolerance,maxiter);
  buckling_beam_solver.solve();
  buckling_beam_solver.write_to_file(filename);

  cx_vec eigval;
  cx_mat eigvec;
  eig_gen(eigval,eigvec,A);
  cout << "\n With Lapack:"<<endl;
  cout << eigval<<endl;
  cout << eigvec << endl;

  return 0;
}

dmat tridiagonal_Toeplitz_matrix(int N, double lower, double center, double upper){
  vec A1 = zeros<vec>(N-1);
  A1(0) = center;
  A1(1) = lower;
  dmat A = toeplitz(A1);
  return A;
}
