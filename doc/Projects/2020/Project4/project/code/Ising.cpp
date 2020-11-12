#include"Ising.hpp"

Ising_2d::Ising_2d(unsigned L,unsigned seed){
  m_L = L;
  rng.seed(seed);
  S = new int*[m_L+2];
  for(int i=0;i<m_L+2;i++){
    S[i] = new int[L+2];
  }
}


Ising_2d::~Ising_2d(){
  for(int i=0;i<m_L+2;i++){delete[] S[i];}
  delete[] S;
}


void Ising_2d::Initialize(double T){
  cout << "Initialize for T = " << T << "\n\t" << flush;
  m_T = T;
  for(int i=-8;i<=8;i+=4){
    w[i+8] = exp(-i/m_T);// precalculate w values
  }
  E = M = 0;
  if(T<1.5){AlignSpins();}
  else{RandomizeSpins();}
//  cout << m_L << endl;

  for(unsigned x=1;x<=m_L;x++){
    for(unsigned y=1;y<=m_L;y++){
//      cout << S(y,x) << endl;
      M += (double) S[y][x];
      E -= (double) S[y][x]*(S[y-1][x] + S[y][x-1]);
    }
  }
//  cout << "\tSpin matrix:\n" << S << flush;
  cout << "\tE  = " << E << " M = " << M << endl;
}


void Ising_2d::Metropolis(){
  int ix,iy;
  for(int i=0;i<m_L*m_L;i++){ // One Monte Carlo cycle
    ix = (int) (uniform(rng)*m_L+1); // +1 because matrix starts at 1
    iy = (int) (uniform(rng)*m_L+1); //   0 is periodic boundary
    int DE = 2*S[iy][ix]*
      (S[iy][ix+1]+
       S[iy][ix-1]+
       S[iy+1][ix]+
       S[iy-1][ix]);
    //if(uniform(rng) <= w(DE+8)){
    bool flip = uniform(rng) <= w[DE+8];
      S[iy][ix] *= -1*(flip*2-1); // flip if flip==1 (true)
      for(int j=1;j<m_L+1;j++){ // Make sure boundaries are correct
        S[0][j]     = S[m_L][j];
        S[m_L+1][j] = S[1][j];
        S[j][0]     = S[j][m_L];
        S[j][m_L+1] = S[j][1];
      }
      M += (double) 2*S[iy][ix]*flip; // zero if flip==0 (false)
      E += (double) DE*flip;
//      cout << E << " " << flush;

  }
}


void Ising_2d::AlignSpins(){
  cout << "AlignSpins" << endl;
  for(int i=0;i<m_L+2;i++){ // Fill S with +1
    for(int j=0;j<m_L+2;j++){
      S[i][j] = +1;
    }
  }
  S[0][0] = S[0][m_L+1] = S[m_L+1][0] = S[m_L+1][m_L+1] = 0;
}


void Ising_2d::RandomizeSpins(){
  cout << "RandomizeSpins" << endl;
  for(int i=0;i<m_L+2;i++){ // Fill S with random +/- 1 
    for(int j=0;j<m_L+2;j++){
      S[i][j] = (uniform(rng)>0.5)*2 - 1;
    }
  }
  for(int j=1;j<m_L+1;j++){ // Set periodic boundaries
    S[0][j]     = S[m_L][j];
    S[m_L+1][j] = S[1][j];
    S[j][0]     = S[j][m_L];
    S[j][m_L+1] = S[j][1];
  }
  S[0][0] = S[0][m_L+1] = S[m_L+1][0] = S[m_L+1][m_L+1] = 0;
}
