#include "EigenValueSolver.hpp"

void Jacobi_rotation::init(dmat A, double tolerance, int maxiter){
  cout << "init" << endl;
  initialize(A);
  m_tolerance = tolerance;
  m_maxiter = maxiter;
}

void Jacobi_rotation::solve(){
  m_iterations = 0;
  m_maxnondiag = 1; // must simply be higher than tolerance
  while(m_maxnondiag > m_tolerance && m_iterations <= m_maxiter){
    max_offdiag();
  //  cout << m_maxnondiag <<endl;
    rotation_angle();
  //  cout << m_s << m_c<<endl;
    rotate();
    m_iterations++;
  //  cout<<m_A<<endl;
  }
  cout << m_iterations<<" iteratons"<<endl;
}

void Jacobi_rotation::max_offdiag(){
  cout << "max_offdiag" << endl;
  m_maxnondiag = 0;
  for(int i=0; i<m_N-1; i++){
    for(int j=i+1; j<m_N-1; j++){
      double max = abs(m_A(i,j));
      if(max > m_maxnondiag){
        m_k = i;
        m_l = j;
        m_maxnondiag = max;
      }
    }
  }
}

void Jacobi_rotation::rotation_angle(){
  cout << "rotation_angle" << endl;
  double tau = (m_A(m_l,m_l)-m_A(m_k,m_k))/(2*m_A(m_k,m_l));
  double t;
  if(tau >= 0){
    t = 1.0/(tau + sqrt(1.0 + tau*tau));
  }
  else{
    t = -1.0/(-tau + sqrt(1.0 + tau*tau));
  }
  m_c = 1.0/sqrt(1 + t*t);
  m_s = t*m_c;
}

void Jacobi_rotation::rotate(){
  cout << "rotate"<<endl;
  double a_kk, a_ll, a_kl, a_ik, a_il, r_ik, r_il;
  a_kk = m_A(m_k,m_k);
  a_ll = m_A(m_l,m_l);
  a_kl = m_A(m_k,m_l);
  m_A(m_k,m_k) = a_kk*m_c*m_c - 2*a_kl*m_c*m_s + a_ll*m_s*m_s;
  m_A(m_l,m_l) = a_ll*m_c*m_c + 2*a_kl*m_c*m_s + a_kk*m_s*m_s;
  m_A(m_k,m_l) = 0.0;
  m_A(m_l,m_k) = 0.0;
  for(int i=0; i<m_N-1; i++){
    if(i != m_k && i != m_l){
      a_ik = m_A(i,m_k);
      a_il = m_A(i,m_l);
      m_A(i,m_k) = a_ik*m_c - a_il*m_s;
      m_A(i,m_l) = a_il*m_c + a_ik*m_s;
      m_A(m_k,i) = m_A(i,m_k);
      m_A(m_l,i) = m_A(i,m_l);
    }
    // New eigenvectors as column vectors
    r_ik = m_R(i,m_k);
    r_il = m_R(i,m_l);
    m_R(i,m_k) = m_c*r_ik - m_s*r_il;
    m_R(i,m_l) = m_c*r_il + m_s*r_ik;
  }
}
