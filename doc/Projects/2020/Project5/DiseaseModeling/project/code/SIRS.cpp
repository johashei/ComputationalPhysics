#include "SIRS.hpp"


/* Constructor: constructs the state of the system at time t=0. Can optionally
set a specific seed for the random number generator.
  Arguments:  int S : Initial number of susceptible
              int I : Initial number of infected
              int R : Initial number of recovered
              int seed (optional) : seed for the random number generator
*/
SIRS::SIRS(int S, int I, int R, int seed/*=1*/){
  cout << "seed = " << seed << endl;
  rng.seed(seed);
  m_state[0] = S;
  m_state[1] = I;
  m_state[2] = R;
  m_state[3] = S+I+R;
}

/* The function calculate_dt calculates the time step such that at most one
individual will move from one group to another. This allowes this number to be
interpreted as a transition probability.
Called by SIRS::Solve_MC
*/
void SIRS::calculate_dt(){
  dt = *min_element(m_invMaxr.begin(),m_invMaxr.end());
  cout << dt << endl;
}

/* The function Solve_MC simulates the system for the given time and writes the
results to the specified file, as tab separated values.
When adding a transition, it should not be necessary to modify this function.
  Arguments:  double final_time : length of the simulation
              int filestep : [NOT IMPLEMENTED] number of timesteps between each write to file
              string outfilename : name of the file to generat
*/
void SIRS::Solve_MC(double final_time, int filestep, string outfilename)
{
  for(auto& x: m_parameters){
    cout << x.first << " : " << x.second << endl;
  }
  calculate_dt();
  int MCs = int(final_time/dt);
  cout << MCs << endl;

  outfile.open(outfilename);
  outfile << m_time << "\t" << m_state[0] << "\t" << m_state[1] << "\t" << m_state[2] << "\n";
  for(int cycle=0;cycle<MCs;cycle++){
    for(int i=0;i<m_functions.size();i++){
      (this->*(m_functions[i]))();
    }
    m_state[3] = m_state[0]+m_state[1]+m_state[2]; // Update N
    m_time += dt;
    outfile << m_time << "\t" << m_state[0] << "\t" << m_state[1] << "\t" << m_state[2] << "\n";
  }
  outfile.close();
}



/*******************************************************************************

The following functions initialize and describe the transitions of the system.
When adding a transition, these two functions should be placed here.
*/


/* The initializer for the basic SIRS model. This should should be called first
after the constructor.
Initialization functions for other transitions should follow this structure.
  Arguments:  double rStoI : transmission parameter
              double rItoR : recovery parameter
              double rRtoS : immunity loss parameter
*/
void SIRS::Basic_SIRS(double rStoI, double rItoR, double rRtoS){
  // Add the parameters to the map m_parameters
  m_parameters["a"]=rStoI;
  m_parameters["b"]=rItoR;
  m_parameters["c"]=rRtoS;
  // Add the address of the function describing the transitions to the vector m_functions
  m_functions.push_back(&SIRS::Basic_SIRS);
  // Add the inverse of the max value of the transition to the vector m_invMaxr
  m_invMaxr.push_back(4/(m_parameters["a"]*m_state[3]));
  m_invMaxr.push_back(1/(m_parameters["b"]*m_state[3]));
  m_invMaxr.push_back(1/(m_parameters["c"]*m_state[3]));
}
/* Function describing the transitions for the basic SIRS model.
Functions describing other transitions should follow this structure.
*/
void SIRS::Basic_SIRS(){
  // Randomly select the transitions according to the transition probabilities
  bool StoI = uniform(rng) < m_parameters["a"]/m_state[3] *m_state[0]*m_state[1]*dt;
  bool ItoR = uniform(rng) < m_parameters["b"]*m_state[1]*dt;
  bool RtoS = uniform(rng) < m_parameters["c"]*m_state[2]*dt;
  // Update the system with the selected transitions
  m_state[0] -= StoI; m_state[1] += StoI;
  m_state[1] -= ItoR; m_state[2] += ItoR;
  m_state[2] -= RtoS; m_state[0] += RtoS;
}



/* The initializer for seasonal variation.
  Arguments:  double amplitude : maximum deviation from the mean transmission rate
              double omega : angular frequency of the oscillation
*/
void SIRS::Seasonal_variation(double amplitude, double omega){
  m_parameters["A"]=amplitude;
  m_parameters["w"]=omega;
  // Replace the Basic_SIRS function by the Seasonal_variation function
  m_functions[0] = (&SIRS::Seasonal_variation);
  m_invMaxr[0] = 4/((m_parameters["a"]+m_parameters["A"])*m_state[3]);
}
/* The transitions for the basic SIRS model, modified with a periodic transmission
parameter.
*/
void SIRS::Seasonal_variation(){
  bool StoI = uniform(rng) < (m_parameters["A"] * cos(m_parameters["w"]*m_time)
                      + m_parameters["a"])/m_state[3] *m_state[0]*m_state[1]*dt;
  bool ItoR = uniform(rng) < m_parameters["b"]*m_state[1]*dt;
  bool RtoS = uniform(rng) < m_parameters["c"]*m_state[2]*dt;
  m_state[0] -= StoI; m_state[1] += StoI;
  m_state[1] -= ItoR; m_state[2] += ItoR;
  m_state[2] -= RtoS; m_state[0] += RtoS;
}



/* The initializer for vital dynamics.
  Arguments:  double rBirth : birth rate
              double rDeath : rate of death not related to the disease
              double rIDeath : rate of death caused by the disease
*/
void SIRS::Vital_dynamics(double rBirth, double rDeath, double rIDeath){
  m_parameters["e"]=rBirth;
  m_parameters["d"]=rDeath;
  m_parameters["dI"]=rIDeath;
  m_functions.push_back(&SIRS::Vital_dynamics);
  m_invMaxr.push_back(1/(m_parameters["e"]*m_state[3]));
  m_invMaxr.push_back(1/(m_parameters["d"]*m_state[3]));
  m_invMaxr.push_back(1/(m_parameters["dI"]*m_state[3]));
}
/* The transitions for vital dynamics. Individuals are added to S by birth, and
removed from all categories by death. Infected individuals can also die from the
disease.
*/
void SIRS::Vital_dynamics(){
  m_state[0] += uniform(rng) < m_parameters["e"]*m_state[3]*dt; // births
  m_state[1] -= uniform(rng) < m_parameters["dI"]*m_state[1]*dt; // deaths from disease
  for(int i=0;i<3;i++){
    m_state[i] -= uniform(rng) < m_parameters["d"]*m_state[i]*dt; // deaths from other causes
  }
}



/* The initializer for vaccination.
  Arguments:  double (*rate)(double,double*) : function describing the vaccination
                rate. Arguments other than time must be passed as a single double pointer.
              double* args : pointer to the arguments other than time to pass to
                the function rate.
*/
void SIRS::Vaccination(double (*rate)(double time, double* args), double* args){
  rVacc = rate;
  m_VaccArgs = args; // Other arguments to pass to the function rate.
  m_functions.push_back(&SIRS::Vaccination);
  // WARNING : Max vaccination rate is not known. If it becomes larger than all
  //         : other transition rates, the algorithm won't be able to keep up.
  //    TODO : Add exception handling for this case.

}
/* The transitions for vacination. Vaccinated individuals move directly from S to R.
*/
void SIRS::Vaccination(){
  // If S=0, there is no one to vaccinate
  bool StoR = uniform(rng) < rVacc(m_time,m_VaccArgs)*(m_state[0]>0)*dt;
  m_state[0] -= StoR;
  m_state[2] += StoR;
}
