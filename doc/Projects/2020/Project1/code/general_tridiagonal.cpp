#include "general_tridiagonal.hpp"
#include <iostream>
#include <cmath>

using namespace std;

void general_tridiagonal::Initialize(int N){
  m_a = new double [N];
  m_b = new double [N];
  for(int i = 0 ; i < m_N ; i++){
    m_a[i] = 0.0;
    m_b[i] = 0.0;
    m_c[i] = 0.0;
  }
}

void general_tridiagonal::solve(){
  double * btilde = new double [N]; // function variable
  double * ftilde = new double [N]; // function variable
  btilde[0] = m_b[0];
  ftilde[0] = m_f[0];

  // Forward substitution
  for(int i=1 ; i<m_N ; i++){
    btilde[i] = m_b[i] - m_a[i]*m_c[i-1]/btilde[i-1];
    ftilde[i] = m_f[i] - m_a[i]*ftilde[i-1]/btilde[i-1];
  }
  //Backward substitution
  m_u[m_N] = ftilde[m_N]/btilde[m_N];
  for(int i=m_N ; i>0 ; i--){
    m_u[i-1] = (ftilde[i-1] - m_c[i-1]*m_u[i])/btilde[i-1]
  }
}
