#ifndef PhysicsSimulator_HPP
#define PhysicsSimulator_HPP
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <string>
#include <cstdio>
#include <ctime>
#include "physics_object.hpp"

using namespace std;

class PhysicsSimulator{

private:
  bool m_CoM; // True : move to CoM system, False : stay in input system.
  unsigned int m_N_objects; // number of moving objects to simulate
  unsigned int m_N_fixed; // number of fixed objects in the system
  vector<PhysicsObject> m_initial; // Vector to store initial conditions
  vector<PhysicsObject> m_fixed; // Vector to store fixed objects. These are not updated during the simulation.
  PhysicsObject *m_objects; // array to fill during simulation. This will contain all positions and velocities.
  double m_simulation_time; // End time of the simulation, start time = 0.
  unsigned int m_N_steps; // number of time steps
  double m_gen_rel; // Correction due to general relativity
  double m_h; // timestep

  vec3 acceleration(unsigned int k, unsigned int i);
  /* Return acceleration of object i at time k*m_h due to gravity from the other
  objects in the system. */

  void velocityVerlet();
  /* Calculate the positions and velocities of all objects using the velocity
  Verlet method. Since the velocity depends on acceleration calculated in the
  same iteration, all positions must be calculated before the velocities.
  Calls PhysicsSimulator::acceleration */

  void Euler();
  /* Calculate the positions and velocities of all objects using the forward
  Euler algorithm.
  Calls PhysicsSimulator::acceleration*/

  void EulerCromer();
  /* Calculate the positions and velocities of all objects using the
  Euler-Cromer algorithm
  Calls PhysicsSimulator::acceleration*/

  void set_initial_conditions();
  /*Find and move to the center of mass frame, and initialize m_objects with
  the position, velocity and mass of each object at time t=0*/

public:
  ~PhysicsSimulator(){cout << "destructor"<<endl;delete [] m_objects;}
  /* Destructor: m_objects allocated in PhysicsSimulator::set_parameters */

  PhysicsSimulator();
  /* Constructor: initialized m_N_objects, m_N_fixed, m_gen_rel and m_CoM */

  void add_fixed_object(PhysicsObject object);
  /* Add a single fixed point particle object to m_fixed. Fixed objects are not
  updated during the simulation. Adding one or more fixed objects also fixes the
  coordinate system. */
  void add_object(PhysicsObject object);
  /* Add single point particle object to m_initial. */

  void add_objects(unsigned int N, const char* posvelFilename, const char* massFilename);
  /*Read position, velocity and mass data for N objects from files and call
  PhysicsSimulator::add_object for each set. The files should be formatted as
  follows:
  posvelFilename: "X Y Z VX VY VZ \n"
  mass: "M \n"
  All additional information should be written after the values. (Sorry)
  */

  void set_parameters(double simulation_time, unsigned int steps, bool gen_rel);
  /*Set the simulation time, the number of time steps to use, and wather to
  include the correction from general relativity. */

  PhysicsObject get_object(unsigned k,unsigned i);
  /* Get a copy of single object i at a single time k*m_h */

  void get_shape(unsigned int shape[2]);
  /* Fill shape with the number of objects and time points */

  void run(string method);
  /* Run the simulation by calling PhysicsSimulator::set_initial_conditions and
  the chosen medtod. */

  void write_to_file(string filename, unsigned int step);
  /*Create a file and write the time and positions and velocities for every
  object with a time point step given by step.*/

  vec3 perihelion(unsigned int k, unsigned int i);
  /* Find the position of the perihelion immediately before the time k*m_h */

};

inline PhysicsObject PhysicsSimulator::get_object(unsigned k, unsigned i){
  return(m_objects[k*m_N_objects + i]);
}

#endif
