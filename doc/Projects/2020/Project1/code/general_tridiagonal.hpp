#ifndef GENERAL_TRIDIAGONAL_HPP
#define GENERAL_TRIDIAGONAL_HPP

class general_tridiagonal {
private:
  int m_N;
  // linear equation is on the form Av = \tilde{b}
  double * m_a;
  double * m_b;
  double * m_c;
  //double * m_btilde;
  double * m_u; // solution
  double * m_f;

public:
  void solve();
}

#endif
