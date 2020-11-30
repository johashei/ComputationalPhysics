#ifndef DEsolver_HPP
#define DEsolver_HPP
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class DEsolver{
private:
protected:
  double m_t;
  vec m_state;
  vec (*coupledDEs)(vec state,double t);
public:
  DEsolver(vec (*f)(vec,double)){m_t = 0; coupledDEs = f;}
  //~DEsolver(){delete [] state;}
  void Set_initial_conditions(vec init_conditions){m_state = init_conditions;}
  vec Get_state(){return m_state;}
  double Get_time(){return m_t;}

};

class RK4 : public DEsolver{
private:
public:
  RK4(vec (*f)(vec,double)):DEsolver(f){;}
  void Advance(double dt); // One time step with the RK4 algorithm
  void Euler(double dt);
};

#endif
