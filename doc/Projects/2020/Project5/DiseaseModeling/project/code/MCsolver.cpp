#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <omp.h>

using namespace std;

int main(int argc, const char* argv[]){
  double a = atof(argv[1]);
  double b = atof(argv[2]);
  double c = atof(argv[3]);
  int S = atoi(argv[4]);
  int I = atoi(argv[5]);
  int R = atoi(argv[6]);
  int MCs = atoi(argv[7]);

  int N = S+I+R;
  double t = 0;
  double dt = min(4/(a*N), min(1/(b*N), 1/(c*N)));

  mt19937_64 rng;
  uniform_real_distribution<double> uniform = std::uniform_real_distribution<double>(0.0,1.0);

  ofstream outfile;
  outfile.open("bA.tsv");

  bool StoI, ItoR, RtoS;
  for(int cycle=0;cycle<MCs;cycle++){
    StoI = uniform(rng) < a*S*I*dt/N;
    ItoR = uniform(rng) < b*I*dt;
    RtoS = uniform(rng) < c*R*dt;
    I += StoI; S -= StoI;
    R += ItoR; I -= ItoR;
    S += RtoS; R -= RtoS;
    t += dt;
    outfile << t << "\t" << S << "\t" << I << "\t" << R << "\n";
  }
  outfile.close();

  return 0;
}
