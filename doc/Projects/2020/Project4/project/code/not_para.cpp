#include <fstream>
#include <iomanip>
#include <ctime>
#include "Ising.hpp"


/*
Unparallelized program for a single temperature. Record E, M, and the number of
accepted configurations and energy after each cycle.
*/

int main(int argc, char const *argv[]){

// Parameters:
  // Side length L
  int L = atoi(argv[1]);
  // Number of Monte Carlo cycles
  int MCS = (int) atof(argv[2]);
  // Define T values
  double T = atof(argv[3]);
  // Initial spin state
  string spin0 = argv[4];
  // Name of file to generate
  string outfilename = argv[5];



  clock_t start_time, run_time;

  int N = L*L; // number of spins;

  ofstream outfile;
  size_t extention_start = outfilename.find_last_of(".");
  outfile.open(outfilename);
  outfile << setw(5) << L << setw(15) << MCS << setw(5) << T << setw(10) << spin0 << "\n";

  // Observables:
  double E,M,N_accepted_configs;

  Ising_2d lattice(L,1);
  lattice.Initialize(T,spin0);

  start_time = clock();
  for(int cycle=0;cycle<MCS;cycle++){
    N_accepted_configs = lattice.Metropolis();
    E = lattice.Get_Energy();
    M = lattice.Get_Magnetization();

    outfile << setw(15) << setprecision(8) << E << "\t";
    outfile << setw(15) << setprecision(8) << M << "\t";
    outfile << setw(15) << setprecision(8) << N_accepted_configs << "\n";
  }
  run_time = (double) (clock() - start_time)/CLOCKS_PER_SEC;
  cout << "Runtime = " << run_time << "s." << endl;

  cout << "Generated file " << outfilename << endl;
  outfile.close();

  return 0;
}
