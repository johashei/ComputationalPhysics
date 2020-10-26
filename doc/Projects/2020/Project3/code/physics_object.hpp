#ifndef PhysicsObject_HPP
#define PhysicsObject_HPP
#include <cmath>
#include "vec3.h"

using namespace std;

class PhysicsObject{ // Fully inlined => header only
private: // The properties of the object
  vec3 pos; // Position of the object
  vec3 vel; // Velocity of the object
  double mass;  // Mass of the object

public: // Setting and getting the properties
  PhysicsObject() {mass = 1;} // Constructor : set pos and vel to 0, mass to 1.
  PhysicsObject(vec3 p, vec3 v, double m); // Constructor with specified properties
  vec3 Get_pos() { return pos; } // Return position as vec3
  double Get_pos(unsigned i) { return pos[i]; } // Return one component of position
  vec3 Get_vel() { return vel; } // Return velocity as vec3
  double Get_vel(unsigned i) { return vel[i]; } // Return one component of velocity
  double Get_mass() { return mass; } // Return mass

  void Set_pos(unsigned i,double value){ pos[i] = value; } // Set one component of position
  void Set_pos(vec3 vector){ pos = vector; } // Set position vector
  void Set_vel(unsigned i,double value){ vel[i] = value; } // Set one component of velocity
  void Set_vel(vec3 vector){ vel = vector; } // Set velocity vector
  void Set_mass(double value){ mass = value; } // Set mass 
};

inline PhysicsObject::PhysicsObject(vec3 p, vec3 v, double m){
  pos = p;
  vel = v;
  mass = m;
}

#endif
