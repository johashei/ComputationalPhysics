#include "physics_simulator.hpp"
#include "tests.hpp"
/*
Use the PhysicsSimulator class to study the Earth-Sun, Earth-Sun-Jupiter and
Mercury-Sun systems.
Difference between Verlet and Euler.
Conservation of energy and angular momentum for different orbits */

int main(int argc, char const *argv[]){

  string outfilename = string(argv[1]);

  vec3 r0Sun(0,0,0);
  vec3 v0Sun(0,0,0);
  double mSun = 1;
  PhysicsObject Sun(r0Sun,v0Sun,mSun);

  vec3 r0Earth(1,0,0);
  vec3 v0Earth(0,2*M_PI,0); // For circular orbit
  double mEarth = 6e24/2e30;
  PhysicsObject Earth(r0Earth,v0Earth,mEarth);

  vec3 r0Mercury(0.3075,0,0);
  vec3 v0Mercury(0,12.44,0);
  double mMercury = 3.3e23/2e30;
  PhysicsObject Mercury(r0Mercury,v0Mercury,mMercury);

  vec3 r0Jupiter(5.20,0,0);
  vec3 v0Jupiter(0,sqrt(4*M_PI*M_PI/5.20),0); // For circular orbit
  double mJupiter = 1.9e27/2e30*100;
  PhysicsObject Jupiter(r0Jupiter,v0Jupiter,mJupiter);


  unsigned int N = 1e8;
  int T = 1;
  PhysicsSimulator TwobodySystem;
  TwobodySystem.add_fixed_object(Sun);
  TwobodySystem.add_object(Mercury);
  //TwobodySystem.add_object(Earth);
  //TwobodySystem.add_object(Jupiter);
  TwobodySystem.set_parameters(T,N,false);
  TwobodySystem.run("velocityVerlet");
  Tests tester(TwobodySystem);
  tester.test_energy_conservation(1e-9);
  tester.test_angular_momentum_conservation(1e-9);
  //vec3 perihelion;// = TwobodySystem.perihelion(3e4,1);
  //cout<< "perihelion (t=0) = " << perihelion <<" = ("<<perihelion.length()<<", "<< atan(perihelion[1]/perihelion[0])*180/M_PI*3600 <<"\")."<<endl;
  vec3 perihelion = TwobodySystem.perihelion(N-1,1);
  cout<< "perihelion = " << perihelion <<" = ("<<perihelion.length()<<", "<< atan(perihelion[1]/perihelion[0])*180/M_PI*3600 <<"\")."<<endl;
  TwobodySystem.write_to_file(outfilename,100);
  cout << "end of program" << endl;
  return 0;
}
