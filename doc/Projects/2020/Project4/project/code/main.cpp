#include <fstream>
#include <iomanip>
#include <omp.h>
#include "Ising.hpp"

#define MAX_THREADS 8

/*
Loop over temperature and number of Monte Carlo cycles. Record expectation values
and number of accepted configurations
*/

int main(int argc, char const *argv[]){

// Parameters:
  // Side length L
  int L = atoi(argv[1]);  //2
  // Define MCS values
  int log10MCS_start = atoi(argv[2]);   //2;
  int log10MCS_end = atoi(argv[3]);     //6;
  double log10MCS_step = atof(argv[4]); //0.5;
  // Define T values
  double T_start = atof(argv[5]); //0.1;
  double T_end = atof(argv[6]);   //3;
  double T_step = atof(argv[7]);  //0.1;
  // Define outfilename
  string outfilename = argv[8]; //"2x2_var_T_MCS_exp.txt"
  // Define initial spin state
  string spin0 = argv[9]; //"default"


  omp_set_num_threads(MAX_THREADS);
  double start_time, run_time;

  int N = L*L; // number of spins;

  int N_MCSvalues = (int) round((log10MCS_end-log10MCS_start)/log10MCS_step) + 1; // include both endpoints
  // round is necessary to avoid error when converting double to int
  int N_Tvalues = (int) round((T_end-T_start)/T_step) + 1; // include both endpoints

  cout << N_MCSvalues << " " << N_Tvalues << endl;

  ofstream outfile,Efile;
  size_t extention_start = outfilename.find_last_of(".");
  string Efilename = outfilename.substr(0,extention_start)
                    +(string)("_E.")
                    +outfilename.substr(extention_start+1);
  outfile.open(outfilename);


  double E[(int)pow(10,log10MCS_end)]; // Cannot init with different number of values for each iteration.

  start_time = omp_get_wtime();
  #pragma omp parallel private(E)
  {
    int id = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    if(id==0){
      printf("num_threads = %d \n",numthreads);
      cout << "Progressbar |" << setw(N_Tvalues*N_MCSvalues) << "|\n" << "            " << flush;
    }
    // Variables:
    double T;
    int MCS;
    // Observables:
    double EE = 0;
    double EM = 0;
    double EE2 = 0;
    double EM2 = 0;
    double EC,EX;
    int N_accepted_configs = 0;

    #pragma omp barrier
    #pragma omp for collapse(2)
    for(int iT=0;iT<N_Tvalues;iT++){
      for(int i=0;i<N_MCSvalues;i++){
        T = T_start + iT*T_step;
        MCS = pow(10,log10MCS_start + i*log10MCS_step);

        Ising_2d lattice(L,id+1);
        lattice.Initialize(T,spin0);

        for(int cycle=0;cycle<MCS;cycle++){
          N_accepted_configs += lattice.Metropolis();

          //bool stable = (cycles > stabilize);
          EE += lattice.Get_Energy();
          EM += abs(lattice.Get_Magnetization());
          EE2 += pow(lattice.Get_Energy(),2);
          EM2 += pow(lattice.Get_Magnetization(),2);
          E[cycle] = lattice.Get_Energy();

        }
        EE = EE/MCS;
        EM = EM/MCS;
        EM2 = EM2/MCS;
        EE2 = EE2/MCS;
        EC = (EE2 - EE*EE)/(T*T);
        EX = (EM2 - EM*EM)/T;

        #pragma omp critical
        {
          outfile << setw(15) << setprecision(8) << T << "\t";
          outfile << setw(15) << setprecision(8) << MCS << "\t";
          outfile << setw(15) << setprecision(8) << N_accepted_configs << "\t";
          outfile << setw(15) << setprecision(8) << EE/N << "\t";
          outfile << setw(15) << setprecision(8) << EM/N << "\t";
          outfile << setw(15) << setprecision(8) << EE2/N << "\t";
          outfile << setw(15) << setprecision(8) << EM2/N << "\t";
          outfile << setw(15) << setprecision(8) << EC/N << "\t";
          outfile << setw(15) << setprecision(8) << EX/N << "\t";
          for(int cycle=0;cycle<MCS;cycle++){
            outfile << E[cycle]/N << " ";
          }
          outfile << "\n";
        }
        cout << "*" << flush; // Update progressbar
      } // Close MCS loop
    } // Close temperature loop
  } // Close parallel
  cout << endl;
  run_time = omp_get_wtime() - start_time;
  cout << "Runtime = " << run_time << "s." << endl;

  cout << "Generated file " << outfilename << endl;
  outfile.close();
  //Efile.close();

  return 0;
}
