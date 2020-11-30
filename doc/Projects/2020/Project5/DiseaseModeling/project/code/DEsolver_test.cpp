#include "DEsolver.hpp"
#include <fstream>

vec DEs(vec SIR, double t);
double a=4, b=1, c=0.5;
int N = 400;

int main(){
  vec Init = {300,100,0};
  RK4 Solver(DEs);
  Solver.Set_initial_conditions(Init);

  ofstream outfile;
  outfile.open("aA.txt");

  for(int i=0;i<50;i++){
    Solver.Advance(0.5);
    outfile << Solver.Get_time() << "\t" << Solver.Get_state().t();
  }
  outfile.close();

  cout << "Expected values at equilibrium :" << endl;
  cout << " S : " << b/a << "\t" << "I : " << (1-b/a)/(1+b/c) << "\t" << "R : " << b/c * (1-b/a)/(1+b/c) << endl;

  return 0;
}

vec DEs(vec SIR, double t){
  vec DSIR(3);
  DSIR(0) = c*SIR(2) - a*SIR(0)*SIR(1)/N;
  DSIR(1) = a*SIR(0)*SIR(1)/N - b*SIR(1);
  DSIR(2) = b*SIR(1) - c*SIR(2);
  return DSIR;
}
