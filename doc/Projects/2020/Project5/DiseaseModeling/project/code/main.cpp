#include "SIRS.hpp"

/* Main program for running the Monte Carlo SIRS model.*/

using namespace std;

double f(double t,double* pf0); // Vaccination

int main(int argc, const char* argv[]){

  double a0 = atof(argv[1]);
  double b = atof(argv[2]);
  double c = atof(argv[3]);
  int S = atoi(argv[4]);
  int I = atoi(argv[5]);
  int R = atoi(argv[6]);
  double Tfinal = atof(argv[7]);
  // Vital dynamics
  double d = atof(argv[8]);
  double dI = atof(argv[9]);
  double e = atof(argv[10]);
  // Seasonal variation
  double A = atof(argv[11]);
  double omega = atof(argv[12]);
  // Vaccination
  double f0 = atof(argv[13]);


  // Generate a different random seed at each run to show seed dependence
  random_device rd;
  int seed = rd();

  SIRS solver(S,I,R,1);
  solver.Basic_SIRS(a0,b,c);
  solver.Vital_dynamics(e,d,dI);
  solver.Seasonal_variation(A,omega);
  solver.Vaccination(f,&f0);
  solver.Solve_MC(Tfinal,1,"SIRSmc.tsv");

  return 0;
}

double f(double t, double* pf0){
  return (*pf0)*(20<t && t<40); // Campaign with constant rate
  //return (*pf0)*(t-20)*(20<t); // Increasing rate
}
