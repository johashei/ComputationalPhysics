#include <fstream>
#include <iomanip>
#include <omp.h>
#include "Ising.hpp"

#define MAX_THREADS 4

/*
Parallelized loop over temperature. Record expectation values
*/

int main(int argc, char const *argv[]){

// Parameters:
  // Side length L
  int L = atoi(argv[1]);
  // Number of Monte Carlo cycles
  int MCS = (int) atof(argv[2]);
  // Discard the first cycles up to
  int equilibration = (int) atof(argv[3]);
  // Define T values
  double T_start = atof(argv[4]);
  double T_end = atof(argv[5]);
  double T_step = atof(argv[6]);
  // Initial spin state
  string spin0 = argv[7];
  // Name of file to generate
  string outfilename = argv[8];



  omp_set_num_threads(MAX_THREADS);
  double start_time, run_time;

  int N = L*L; // number of spins;
  int keptMCS = MCS - equilibration; // Number of cycles used in the calculation of expectation values

  // round is necessary to avoid error when converting double to int
  int N_Tvalues = (int) round((T_end-T_start)/T_step) + 1; // include both endpoints

  ofstream outfile;
  size_t extention_start = outfilename.find_last_of(".");
  outfile.open(outfilename);
  outfile << setw(5) << L << setw(15) << MCS << setw(15) << equilibration << setw(10) << spin0 << "\n";

  start_time = omp_get_wtime();
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    if(id==0){
      printf("num_threads = %d \n",numthreads);
      cout << "Progressbar |" << setw(N_Tvalues) << "|\n" << "            " << flush;
    }
    // Variables:
    double T;
    // Observables:
    double EE = 0;
    double EM = 0;
    double EE2 = 0;
    double EM2 = 0;
    double EC,EX;

    #pragma omp barrier
    #pragma omp for
    for(int iT=0;iT<N_Tvalues;iT++){
        T = T_start + iT*T_step;

        Ising_2d lattice(L,id+1);
        lattice.Initialize(T,spin0);

        for(int cycle=0;cycle<MCS;cycle++){
          lattice.Metropolis();

          bool keepvalue = (cycle > equilibration);
          EE += lattice.Get_Energy()*keepvalue;
          EM += abs(lattice.Get_Magnetization())*keepvalue;
          EE2 += pow(lattice.Get_Energy(),2)*keepvalue;
          EM2 += pow(lattice.Get_Magnetization(),2)*keepvalue;
        }
        EE = EE/keptMCS;
        EM = EM/keptMCS;
        EM2 = EM2/keptMCS;
        EE2 = EE2/keptMCS;
        EC = (EE2 - EE*EE)/(T*T);
        EX = (EM2 - EM*EM)/T;

        #pragma omp critical
        {
          outfile << setw(15) << setprecision(8) << T << "\t";
          outfile << setw(15) << setprecision(8) << EE/N << "\t";
          outfile << setw(15) << setprecision(8) << EM/N << "\t";
          outfile << setw(15) << setprecision(8) << EE2/N << "\t";
          outfile << setw(15) << setprecision(8) << EM2/N << "\t";
          outfile << setw(15) << setprecision(8) << EC/N << "\t";
          outfile << setw(15) << setprecision(8) << EX/N << "\n";
        }
        cout << "*" << flush; // Update progressbar

    } // Close temperature loop
  } // Close parallel
  cout << endl;
  run_time = omp_get_wtime() - start_time;
  cout << "Runtime = " << run_time << "s." << endl;

  cout << "Generated file " << outfilename << endl;
  outfile.close();

  return 0;
}
