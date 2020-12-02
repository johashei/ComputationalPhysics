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
  double d = atof(argv[8]);
  double dI = atof(argv[9]);
  double e = atof(argv[10]);

  int N = S+I+R;
  double t = 0;

  double minTs[] = {4/(a*N), 1/(b*N), 1/(c*N), 1/(e*N), 1/(d*N)};
  double dt = *min_element(minTs,minTs+5);

  mt19937_64 rng;
  uniform_real_distribution<double> uniform = std::uniform_real_distribution<double>(0.0,1.0);

  ofstream outfile;
  outfile.open("bA.tsv");

  bool StoI, ItoR, RtoS, birth, Sdeath, Ideath, Rdeath;
  for(int cycle=0;cycle<MCs;cycle++){
    StoI = uniform(rng) < a*S*I*dt/N;
    ItoR = uniform(rng) < b*I*dt;
    RtoS = uniform(rng) < c*R*dt;
    birth = uniform(rng) < e*N*dt;
    Sdeath = uniform(rng) < d*S*dt;
    Ideath = uniform(rng) < (d+dI)*I*dt;
    Rdeath = uniform(rng) < d*R*dt;
    I += StoI; S -= StoI;
    R += ItoR; I -= ItoR;
    S += RtoS; R -= RtoS;
    S += birth;
    S -= Sdeath; I -= Ideath; R -= Rdeath;
    t += dt;
    outfile << t << "\t" << S << "\t" << I << "\t" << R << "\n";
  }
  outfile.close();

  return 0;
}
