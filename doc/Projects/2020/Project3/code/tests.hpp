#ifndef Tests_HPP
#define Tests_HPP
#include "physics_simulator.hpp"

using namespace std;

class Tests{
private:
  unsigned int m_shape[2]; // array to stor number of time points and objects
  //vector<double> m_distances;
  PhysicsSimulator *m_system; // pointer to the instance of PhysicsSimulator to test

  double kinetic_E(unsigned k);
  /* Return the kinetic energy of m_system at time point k */

  double potential_E_dynamic(unsigned k);
  /* Return the potential energy from the moving objects in m_system at time point k */

  double potential_E_fixed(unsigned k);
  /* Return the potential energy from the fixed objects in m_system at time point k */

  vec3 angular_momentum(unsigned k);
  /* Return the total angular momentum of m_system at time point k */

public:
  Tests(PhysicsSimulator& system);
  /* Constructor: take a reference to PhysicsSimulator instance to initialize *m_system */

  int test_energy_conservation(double tolerance);
  /* Calculate the total energy of m_system at the first and last time points
  and compare. Return 1 if the difference is whithin tolerance, 0 else */

  int test_angular_momentum_conservation(double tolerance);
  /* Calculate the total angular momentum of m_system at the first and last time
  points and compare. Return 1 if the difference is whithin tolerance, 0 else */
};

#endif
