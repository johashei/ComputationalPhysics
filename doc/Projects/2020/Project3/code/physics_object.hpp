#ifndef PhysicsObject_HPP
#define PhysicsObject_HPP
#include <cmath>
#include "vec3.h"

using namespace std;

class PhysicsObject{ // Fully inlined => header only
private:
  vec3 pos;
  vec3 vel;
  double mass;
public:
  PhysicsObject() {mass = 1;} // zero mass not supported
  PhysicsObject(vec3 p, vec3 v, double m);
  vec3 Get_pos() { return pos; } // Return position as vec3
  double Get_pos(unsigned i) { return pos[i]; } // return position component i
  vec3 Get_vel() { return vel; }
  double Get_vel(unsigned i) { return vel[i]; }
  double Get_mass() { return mass; }

  void Set_pos(unsigned i,double value){ pos[i] = value; }
  void Set_pos(vec3 vector){ pos = vector; }
  void Set_vel(unsigned i,double value){ vel[i] = value; }
  void Set_vel(vec3 vector){ vel = vector; }
  void Set_mass(double value){ mass = value; }
};

inline PhysicsObject::PhysicsObject(vec3 p, vec3 v, double m){
  pos = p;
  vel = v;
  mass = m;
}

#endif
