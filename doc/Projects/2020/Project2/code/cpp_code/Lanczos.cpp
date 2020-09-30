#include "EigenValueSolver.hpp"

void Lanczos::init(dmat A, uword M){
  initialize(A);
  m_M = M;
}

void Lanczos::solve(){
  iterate();
  diagonalize(m_M);
}

void Lanczos::iterate(){
  cout << "iterating" << endl;
  m_alpha.zeros(m_M+1);
  m_beta.ones(m_M+1);
  vec r(m_N-1,fill::zeros);r(0) = 1;
  vec q(m_N-1,fill::zeros);
  vec q_old(m_N-1); // q_{k-1}
  vec Aq(m_N-1);
  uword k = 0;

  while(k < m_M){
    q_old = q;
    q = r/m_beta(k);
    Aq = m_A * q;
    k = k+1;
    m_alpha(k) = dot(q,Aq);
    r = Aq - m_alpha(k)*q - m_beta(k-1)*q_old;
    m_beta(k) = norm(r,2);

    //if(k%10==0){ // Print eigenvalues inside loop
    //  diagonalize(k);
    //  cout << "k = " << k << endl;
    //  cout << m_eigvals.subvec(0,4) << endl;
    //}
  }
}

void Lanczos::diagonalize(uword k){
  cout << "Diagonalizing" << endl;
  mat T = zeros<mat>(k,k);
  T.diag(0) = m_alpha.subvec(1,k);
  T.diag(1) = m_beta.subvec(1,k-1);
  T.diag(-1) = m_beta.subvec(1,k-1);
  m_eigvals = eig_sym(T);
}

vec Lanczos::Get_eigvals(uword from, uword to){
  return(m_eigvals.subvec(from,to));
}
