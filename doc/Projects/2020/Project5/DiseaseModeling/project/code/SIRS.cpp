#include "SIRS.hpp"

void SIRS::Solve_MC(double final_time, int filestep, string outfilename){
  int N = m_state[0] + m_state[1] + m_state[2];
  double minTs[] = {
    4/((m_parameters["a/N"]+m_parameters["A"])*N*N),
    1/(m_parameters["b"]*N),
    1/(m_parameters["c"]*N)
  };
  for(auto& x: m_parameters){
    cout << x.first << " : " << x.second << endl;
  }
  dt = *min_element(minTs,minTs+3);
  cout << dt << endl;
  int MCs = int(final_time/dt);
  cout << MCs << endl;

  outfile.open(outfilename);
  for(int cycle=0;cycle<MCs;cycle++){
    SeasonalSIRS();
    Vital_dynamics();
    //Vaccination();
    m_time += dt;
    outfile << m_time << "\t" << m_state[0] << "\t" << m_state[1] << "\t" << m_state[2] << "\n";
  }
  outfile.close();
}

void SIRS::SeasonalSIRS(){
  /* The basic SIRS model, but with periodic transmission parameter. */
  bool StoI = uniform(rng) < (m_parameters["A"] * cos(m_parameters["w"]*m_time) + m_parameters["a/N"])*m_state[0]*m_state[1]*dt;
  bool ItoR = uniform(rng) < m_parameters["b"]*m_state[1]*dt;
  bool RtoS = uniform(rng) < m_parameters["c"]*m_state[2]*dt;
  m_state[0] -= StoI; m_state[1] += StoI;
  m_state[1] -= ItoR; m_state[2] += ItoR;
  m_state[2] -= RtoS; m_state[0] += RtoS;
}

void SIRS::Vital_dynamics(){
  m_state[1] += uniform(rng) < m_parameters["e"]*dt; // births
  m_state[2] -= uniform(rng) < m_parameters["dI"]*m_state[1]*dt; // deaths from disease
  for(int i=1;i<4;i++){
    m_state[i] -= m_parameters["d"]*m_state[i]*dt; // deaths from other causes
  }
}

void SIRS::Vaccination(){
  cout << "vacc" << endl;
  bool StoR = uniform(rng) < rVacc(m_time)*dt;
  cout << StoR << endl;
  m_state[0] -= StoR;
  m_state[2] += StoR;
}
