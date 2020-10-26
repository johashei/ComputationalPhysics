#include "physics_simulator.hpp"
#include "tests.hpp"

/* Simulate the solar system with 8 planets and Pluto */

int main(int argc, char const *argv[]){

  string outfilename = string(argv[1]);

  vec3 r0Sun(0,0,0);
  vec3 v0Sun(0,0,0);
  double mSun = 1;

  const char* pvfile = "posvel.dat";
  const char* mfile = "masses.dat";

  unsigned N = 1e5;
  unsigned T = 10;
  PhysicsSimulator SolarSystem;
  PhysicsObject Sun(r0Sun, v0Sun, mSun);

cout << "add sun"<< endl;
  SolarSystem.add_object(Sun);
cout<<"add planets"<<endl;
  SolarSystem.add_objects(9,pvfile,mfile);
cout<<"parameters"<<endl;
  SolarSystem.set_parameters(T,N,0);
cout<<"run"<<endl;
  SolarSystem.run("velocityVerlet");
cout<<"test"<<endl;
  Tests tester(SolarSystem);
  tester.test_energy_conservation(1e-7);
  tester.test_angular_momentum_conservation(1e-7);
  SolarSystem.write_to_file(outfilename,N/(T*50));

  return 0;
}
