#ifndef EigenValueSolver_HPP
#define EigenValueSolver_HPP
#include <fstream>
#include <armadillo>
#include <iostream>
#include <ctime>

using namespace std;
using namespace arma;

class EigenValueSolver{
private:
protected:
  int m_N;
  ofstream m_ofile;
  dmat m_A; // matrix to diagonalize

public:
  void initialize(dmat A);
  void write_to_file(string filename);

  // These are public to be used by tests. todo : find a better way
  dmat m_R; // matrix for storing eigenvectors
};

class Jacobi_rotation : public EigenValueSolver{
private:
  double m_tolerance;
  int m_iterations, m_maxiter;
  double m_c, m_s;
  void rotation_angle();
  void rotate();
public:
  void init(dmat A, double tolerance, int maxiter);
  void solve();
  vec Get_eigvals(uword from, uword to);


  // These are public to be used by tests. todo : find a better way
  void max_offdiag();
  int m_k, m_l;
  double  m_maxnondiag;
};

class Lanczos : public EigenValueSolver{
private:
  uword m_M;
  vec m_r;
  vec m_alpha, m_beta;
  vec m_eigvals;
  void iterate();
  void diagonalize(uword k);
public:
  void init(dmat A, uword M);
  void solve();
  vec Get_eigvals(uword from, uword to);
};



#endif
