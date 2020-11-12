#include <fstream>
#include <iomanip>
#include <omp.h>
#include "Ising.hpp"

#define MAX_THREADS 8

int main(int argc, char const *argv[]){
  /*double T,T_start,T_end,T_step,Tc,E,M;
  unsigned MCS,L,N_Tvalues;
  Tc = 2.609;
  MCS = 1e5;
  L = atoi(argv[1]); // side length of the spin lattice.
  T_start = atoi(argv[2]); // T in units of kT/J
  T_end = atoi(argv[3]);
  T_step = atoi(argv[4]);

  N_Tvalues = (int) (T_end-T_start)/T_step + 1;
  vec EE(N_Tvalues);
  vec EM(N_Tvalues);*/
  unsigned L = atoi(argv[1]);
  unsigned MCS = atof(argv[2]);
  //double T;// = atoi(argv[3]);

  double T_start = 0.1;
  double T_end = 3;
  double T_step = 0.1;
  int N_Tvalues = (int) ((T_end-T_start)/T_step) + 1; // include both endpoints
  cout << N_Tvalues << endl;


  omp_set_num_threads(MAX_THREADS);

  ofstream outfile;
  char outfilename[30];
  sprintf(outfilename,"%dx%d_exp.txt",L,L);
  outfile.open(outfilename);

  //int counter = 0; // There is probably a better way to do this
  //#pragma omp parallel
  //int id = omp_get_thread_num();
  //int numthreads = omp_get_num_threads();

  //double T_i = T_start + T_step*id;

  double avgfactor = (L*L*MCS); // Average per spin;

  #pragma omp parallel
  {
    int id = omp_get_thread_num();
    int numthreads = omp_get_num_threads();
    if(id==0){printf(" num_threads = %d \n",numthreads);}

    double EE = 0;
    double EM = 0;
    double EE2 = 0;
    double EM2 = 0;
    double EC,EX,T;

    #pragma omp for
    for(int iT=0;iT<N_Tvalues;iT++){
      T = T_start + iT*T_step;
      Ising_2d lattice(L,1240);
  //   #pragma omp parallel
  //  {
  //    int id = omp_get_thread_num();
  //    int numthreads = omp_get_num_threads();
  //    if(id==0){printf(" num_threads = %d \n",numthreads);}

      lattice.Initialize(T);
  //    #pragma omp for reduction(+:EE) reduction(+:EM) reduction(+:EE2) reduction(+:EM2)
      for(int cycles=1;cycles<=MCS;cycles++){
        lattice.Metropolis();
        //outfile << setw(15) << setprecision(8) << lattice.Get_Energy()/(L*L) << "\t";
        //outfile << setw(15) << setprecision(8) << lattice.Get_Magnetization()/(L*L) << "\n";
        EE += lattice.Get_Energy();
        EM += abs(lattice.Get_Magnetization());
        EE2 += pow(lattice.Get_Energy(),2);
        EM2 += pow(lattice.Get_Magnetization(),2);
      }
    //}
      EE = EE/avgfactor;
      EM = EM/avgfactor;
      EM2 = EM2/avgfactor;
      EE2 = EE2/avgfactor;
      EC = (EE2 - EE*EE)/(T*T);
      EX = (EM2 - EM*EM)/T;
      #pragma omp critical
      {
        outfile << setw(15) << setprecision(8) << T << "\t";
        outfile << setw(15) << setprecision(8) << EE << "\t";
        outfile << setw(15) << setprecision(8) << EM << "\t";
        outfile << setw(15) << setprecision(8) << EE2 << "\t";
        outfile << setw(15) << setprecision(8) << EM2 << "\t";
        outfile << setw(15) << setprecision(8) << EC << "\t";
        outfile << setw(15) << setprecision(8) << EX << "\n";
      }
    }
  } // Close parallel

  //cout << "<E>/MCS : "<< EE << "\n<M>/MCS : " << EM << endl;
  cout << "Generated file " << outfilename << endl;
  outfile.close();

  return 0;
}
