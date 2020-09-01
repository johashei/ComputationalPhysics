#ifndef GENERAL_TRIDIAGONAL_HPP
#define GENERAL_TRIDIAGONAL_HPP

#include <string>
using namespace std;

class general_tridiagonal {
private:
  int m_N;
  // linear equation is on the form Av = \tilde{b}
  double * m_a;
  double * m_b;
  double * m_c;
  //double * m_btilde;
  double * m_u; // solution
  double * m_g; // function g = h*f values
  

public:
  void Initialize(int N, double* a, double* b, double* c, double f(double x));
  void Solve();
  void Test_analytical();
  void Writefile(string outfilename);

  ~general_tridiagonal(){ // Destructor
    delete [] m_a; delete [] m_b; delete [] m_c; delete [] m_u; delete [] m_g;
  }
};

#endif
