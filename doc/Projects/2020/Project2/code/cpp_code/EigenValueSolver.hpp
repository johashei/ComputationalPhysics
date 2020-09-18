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
  dmat m_A; // matrix to diagonalize
  dmat m_R; // matrix for storing eigenvectors
  int m_N;
  ofstream m_ofile;

public:
  void initialize(dmat A);
  void write_to_file(string filename);
};

class Jacobi_rotation : public EigenValueSolver {
private:
  double m_tolerance, m_maxnondiag;
  int m_iterations, m_maxiter;
  int m_k, m_l;
  double m_c, m_s;
  void max_offdiag();
  void rotation_angle();
  void rotate();
public:
  void init(dmat A, double tolerance, int maxiter);
  void solve();
};


class ArmadilloSolver : public EigenValueSolver {
private:
};

class Test : public EigenValueSolver {
};

#endif
