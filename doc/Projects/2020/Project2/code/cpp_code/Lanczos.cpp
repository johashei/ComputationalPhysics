#include "EigenValueSolver.hpp"

void Lanczos::init(dmat A, uword M){
  initialize(A);
  m_M = M;
  //m_r = r_0;
}

void Lanczos::solve(){
  iterate();
  diagonalize();
}

void Lanczos::iterate(){
  cout << "iterating" << endl;
  m_alpha.zeros(m_M+1);
  m_beta.ones(m_M+1);
  vec r(m_N-1,fill::ones);
  vec q(m_N-1,fill::zeros);
  vec q_old(m_N-1); // q_{k-1}
  vec Aq(m_N-1);
  uword k = 0;

  //cout << m_A << endl;
  //double t; vec w=r; vec v=m_A*w; m_alpha(1)=dot(w,v); v=v-m_alpha(1)*w; m_beta(1)=norm(v,2);k=1;
  while(k < m_M){
    //for(uword i=0;i<m_N-1;i++){
    //  t = w(i);
    //  w(i) = v(i)/m_beta(k);
    //  v(i) = -m_beta(k)*t;
    //}
    //v = v + m_A*w;
    //k = k+1;
    //m_alpha(k) = dot(w,v);
    //v = v - m_alpha(k)*w;
    //m_beta(k) = norm(v,2);
    q_old = q;
    q = r/m_beta(k);
    Aq = m_A * q;
    k = k+1;
    m_alpha(k) = dot(q,Aq);
    r = Aq - m_alpha(k)*q - m_beta(k-1)*q_old;
    m_beta(k) = norm(r,2);
    if(k%10 == 0){
      diagonalize();
      cout << "k = " << k << endl;
      cout << m_eigvals.subvec(0,4) << endl;
    }
    //cin.get();
  }
}

void Lanczos::diagonalize(){
  cout << "Diagonalizing" << endl;
  mat T = zeros<mat>(m_M,m_M);
  T.diag(0) = m_alpha.subvec(1,m_M);
  T.diag(1) = m_beta.subvec(1,m_M-1);
  T.diag(-1) = m_beta.subvec(1,m_M-1);
  m_eigvals = eig_sym(T);
}

vec Lanczos::Get_eigvals(uword from, uword to){
  return(m_eigvals.subvec(from,to));
}
