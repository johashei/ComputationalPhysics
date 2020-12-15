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

  double dt;
  double m_time;
  vector<int> m_state = std::vector<int>(3); // State of the system, as {S,I,R}
  vector<void (SIRS::*)()> m_functions;
  // Functions to call when updating the system:
  void SeasonalSIRS();
  void Vital_dynamics();
  void Vaccination();
    double (*rVacc)(double time, double* args);
    double* m_VaccArgs;
public:
  unordered_map<string,double> m_parameters;
  // Variable names starting with r + [capital letter] represent rates
  SIRS(double rStoI, double rItoR, double rRtoS);
  void Set_state(int S, int I, int R){m_state[0] = S; m_state[1]=I; m_state[2]=R;}
  void Solve_MC(double final_time, int filestep, string outfilename);
  vector<int> Get_state(){return m_state;}
  double Get_time(){return m_time;}
  // Improvements on the basic SIRS modeL: parameter inputs
  void Seasonal_variation(double amplitude, double omega);
  void Vital_dynamics(double rBirth, double rDeath, double rIDeath);
  void Vaccination(double (*rate)(double time, double* args), double* args);
};

inline SIRS::SIRS(double rStoI, double rItoR, double rRtoS){
  m_parameters["a/N"]=rStoI;
  m_parameters["b"]=rItoR;
  m_parameters["c"]=rRtoS;
  m_functions.push_back(&SIRS::SeasonalSIRS);
}

inline void SIRS::Vital_dynamics(double rBirth, double rDeath, double rIDeath){
  m_parameters["e"]=rBirth;
  m_parameters["d"]=rDeath;
  m_parameters["dI"]=rIDeath;
  m_functions.push_back(&SIRS::Vital_dynamics);

}

inline void SIRS::Seasonal_variation(double amplitude, double omega){
  m_parameters["A"]=amplitude;
  m_parameters["w"]=omega;
  // Does not require calling an additional function
}

inline void SIRS::Vaccination(double (*rate)(double time, double* args), double* args){
  rVacc = rate;
  m_VaccArgs = args; // Other arguments to pass to the function rate.
  m_functions.push_back(&SIRS::Vaccination);
}
