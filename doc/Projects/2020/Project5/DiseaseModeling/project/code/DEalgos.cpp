#include "DEsolver.hpp"

/* Update the state of the system with fourth order Runge-Kutta algorithm */
void RK4::Advance(double dt){
  double t15 = m_t + 0.5*dt;
  vec k1 = dt * coupledDEs(m_state,m_t);
  vec k2 = dt * coupledDEs(m_state + 0.5*k1, t15);
  vec k3 = dt * coupledDEs(m_state + 0.5*k2, t15);
  vec k4 = dt * coupledDEs(m_state + k3, m_t+dt);
  m_state = m_state + (k1 + 2*k2 + 2*k3 + k4)/6.0;
  m_t = m_t + dt;
}

void RK4::Euler(double dt){//Working name
  m_state = m_state + coupledDEs(m_state,m_t)*dt;
  m_t = m_t + dt;
}
