#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>
/* Storing parameters in a map allows for more readable code, and ability to
include or not any number of effects*/

using namespace std;

class SIRS{

private:
  mt19937_64 rng;
  uniform_real_distribution<double> uniform = std::uniform_real_distribution<double>(0.0,1.0);
  ofstream outfile;
  void calculate_dt();
  double dt;
  double m_time = 0;
  vector<void (SIRS::*)()> m_functions;
  vector<int> m_state = std::vector<int>(4); // State of the system, as {S,I,R,N}
  vector<double> m_invMaxr;

  /* Functions to call when updating the system. When adding a new transition,
  the function describing it should be declared here. */
  void Basic_SIRS();
  void Seasonal_variation();
  void Vital_dynamics();
  void Vaccination();
    double (*rVacc)(double time, double* args);
    double* m_VaccArgs;

public:
  unordered_map<string,double> m_parameters;
  // Variable names starting with r + [capital letter] represent rates
  SIRS(int S, int I, int R, int seed=1 );
  void Solve_MC(double final_time, int filestep, string outfilename);
  vector<int> Get_state(){return m_state;}
  double Get_time(){return m_time;}

  /* Parameter inputs beyond the basic SIRS model. When adding a new transition,
  its initializer function should be declared here. */
  void Basic_SIRS(double rStoI, double rItoR, double rRtoS);
  void Seasonal_variation(double amplitude, double omega);
  void Vital_dynamics(double rBirth, double rDeath, double rIDeath);
  void Vaccination(double (*rate)(double time, double* args), double* args);
};
