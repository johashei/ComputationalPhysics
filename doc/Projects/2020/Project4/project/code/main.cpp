#include <fstream>
#include <iomanip>
#include <omp.h>
#include "Ising.hpp"

#define MAX_THREADS 8

int main(int argc, char const *argv[]){

  /*if(argc < 5){
    cout << "Bad usage: requires arguments L MCS T_start T_end T_step" << endl;
    exit(1);
  }*/

  omp_set_num_threads(MAX_THREADS);
  double start_time, run_time;

  unsigned L = atoi(argv[1]);
  ///*
  unsigned MCS = atof(argv[2]);
  double T_start = atof(argv[3]);//0.1;
  double T_end = atof(argv[4]);//3;
  double T_step = atof(argv[5]);//0.1;
  int N_Tvalues = (int) ((T_end-T_start)/T_step) + 1; // include both endpoints
  //*/


  //unsigned stabilize = 0.1*MCS;
  //cout << "MCS cutoff " << stabilize << endl;


  ofstream outfile;
  char outfilename[30];
  sprintf(outfilename,"%dx%d_exp.txt",L,L);
  //sprintf(outfilename,"%dx%d_mcs.txt",L,L);
  outfile.open(outfilename);

  start_time = omp_get_wtime();
  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    if(id==0){printf(" num_threads = %d \n",numthreads);}

    double EE = 0;
    double EM = 0;
    double EE2 = 0;
    double EM2 = 0;
    double EC,EX,T=1;

    //unsigned MCS;
    double N = L*L; // number of spins;

    #pragma omp for
    for(int iT=0;iT<N_Tvalues;iT++){
      T = T_start + iT*T_step;
    /*for(int i=2;i<7;i++){
      MCS = pow(10,i);
      avgfactor = (L*L*MCS);
      cout << "MCS = " << MCS << " for thread " << omp_get_thread_num() << endl;*/

      Ising_2d lattice(L,1240);
      lattice.Initialize(T);
      for(int cycles=1;cycles<=MCS;cycles++){
        lattice.Metropolis();

        //bool stable = 1;//(cycles > stabilize); // Don't use the first 10 % of cycles
        EE += lattice.Get_Energy();
        EM += abs(lattice.Get_Magnetization());
        EE2 += pow(lattice.Get_Energy(),2);
        EM2 += pow(lattice.Get_Magnetization(),2);
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
        outfile << setw(15) << setprecision(8) << EE/N << "\t";
        outfile << setw(15) << setprecision(8) << EM/N << "\t";
        outfile << setw(15) << setprecision(8) << EE2/N << "\t";
        outfile << setw(15) << setprecision(8) << EM2/N << "\t";
        outfile << setw(15) << setprecision(8) << EC/N << "\t";
        outfile << setw(15) << setprecision(8) << EX/N << "\n";
      }
    } // Close temperature loop
  } // Close parallel
  run_time = omp_get_wtime() - start_time;
  cout << "Runtime = " << run_time << "s." << endl;

  cout << "Generated file " << outfilename << endl;
  outfile.close();

  return 0;
}
