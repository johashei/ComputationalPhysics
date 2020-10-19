#ifndef PhysicsSimulator_HPP
#define PhysicsSimulator_HPP
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <string>
#include "physics_object.hpp"

using namespace std;

class PhysicsSimulator{

private:
public: // For debugging
  unsigned int m_N_objects; // number of objects to simulate
  vector<PhysicsObject> m_initial; // Vectors to store initial conditions : can change length
  PhysicsObject *m_objects; // arrays to fill during simulation
  double m_simulation_time;
  unsigned int m_N_steps; // number of time steps
  bool m_gen_rel;
  double m_h; // timestep

  //double gravity(PhysicsObject,PhysicsObject); // Gravitational force between two point particle objects
  vec3 acceleration(unsigned k, unsigned i); // Acceleration for object (k,i)

  void velocityVerlet();
  void set_initial_conditions();
  //void scale();
//public:
  PhysicsSimulator();
  void add_object(PhysicsObject object); // Add single point particle object
  void set_parameters(double simulation_time, unsigned int steps, bool gen_rel);
  void run();
  void write_to_file(string filename);
};


#endif
