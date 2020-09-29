#include "EigenValueSolver.hpp"

void EigenValueSolver::initialize(dmat A){
  m_N = A.n_cols + 1; // A is an (n-1)x(n-1)-matrix
  if(m_N-1 != A.n_rows){
    cout << "Bad usage: A should be a square matrix but has size " << size(A) << endl;
  }
  m_A = A;
  m_R.eye(m_N-1,m_N-1);
}


void EigenValueSolver::write_to_file(string filename){
  // Todo: modify to work with Lanczos output 
  cout <<"writing to file "<< filename << endl;// \nA\n"<<m_A<<"\nR\n"<<m_R<<endl;
  m_ofile.open(filename);
  m_ofile << "eigenvalues   &   eigenvectors \n";
  for(int i=0; i<m_N-1; i++){

    //vec r_i = m_R.row(i);
    m_ofile << m_A(i,i) << "   &   " << m_R.col(i) << "\n";
  }
  m_ofile.close();
}
