#include "physics_simulator.hpp"

using namespace std;

int main(){

  vec3 r0Sun(0,0,0);
  vec3 v0Sun(0,0,0);
  double mSun = 1;

  vec3 r0Earth(9.445675412475125E-01,3.223034192685008E-01, -1.992059498358081E-05);
  vec3 v0Earth(-5.831624954084838E-03,1.622307605140163E-02,-1.112702817686854E-06);
  v0Earth = v0Earth*365.25; // Values given in AU/day
  //vec3 r0Earth(1,0,0);
  //vec3 v0Earth(0,2*M_PI,0);
  double mEarth = 6e24/2e30;

  vec3 r0Jupiter(2.536919841360543,-4.451481277859586,-3.827083809378690E-02);
  vec3 v0Jupiter(6.472015893710071E-03,4.096419901354590E-03,-1.618252353577887E-04);
  v0Jupiter = v0Jupiter*365.25;
  double mJupiter = 1.9e27/2e30;


  PhysicsSimulator SolarSystem;
  PhysicsObject Sun(r0Sun, v0Sun, mSun);
  PhysicsObject Earth(r0Earth, v0Earth, mEarth);
  PhysicsObject Jupiter(r0Jupiter, v0Jupiter, mEarth);

  SolarSystem.add_object(Earth);
  SolarSystem.add_object(Sun);
  SolarSystem.add_object(Jupiter);
  SolarSystem.set_parameters(10,500,0);
  SolarSystem.set_initial_conditions();
  cout << "velocityVerlet" << endl;
  SolarSystem.velocityVerlet();
  SolarSystem.write_to_file("classdata.dat");

  return 0;
}
