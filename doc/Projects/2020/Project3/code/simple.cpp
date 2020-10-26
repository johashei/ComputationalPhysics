#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include "vec3.h"

using namespace std;

/* Simple program for a two-body system without object orientation. */

int main(){

  // Set initial conditions:
  vec3 r0Sun(0,0,0);
  vec3 v0Sun(0,0,0);
  double mSun = 1;

  vec3 r0Earth(0.3075,0,0);
  vec3 v0Earth(0,12.44,0);
  cout << r0Earth.dot(v0Earth);

  double mEarth = 3.3e23/2e30;

  // Relativistic correction
  vec3 l = r0Earth.cross(v0Earth);
  vec3 r = r0Earth-r0Sun;
  double c = 63198;
  double gen_rel = 3*l.lengthSquared() / (r.lengthSquared()*c*c);
  cout << gen_rel << endl;

  // Set simulation parameters
  const unsigned N = 1e8;
  double T = 1; // sim time in years
  double h = T/(N-1);
  double hh = h*h;

  vec3 rEarth[3];
  vec3 rSun[3];
  vec3 vEarth[3];
  vec3 vSun[3];
  vec3 rCoM = (mSun*r0Sun + mEarth*r0Earth)/(mSun+mEarth);
  vec3 vCoM = (mSun*v0Sun + mEarth*v0Earth)/(mSun+mEarth);
  cout << rCoM << " " << vCoM << endl;

  rEarth[1] = r0Earth;// - rCoM;
  rSun[1] = r0Sun;// - rCoM;
  vEarth[1] = v0Earth;// - vCoM;
  vSun[1] = v0Sun;// - vCoM;

  // Velocity Verlet
  double GMoE = 4*M_PI*M_PI * mEarth; // ≈ G*M_sun*M_earth
  double mEi = 1/mEarth;
  double mSi = 1/mSun;

  vec3 diff = rEarth[1]-rSun[1];
  double R = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
  vec3 FES = -GMoE*diff/pow(R,3) * (1+gen_rel);
  vec3 FESnew;


  ofstream outfile;
  outfile.open("simple.dat");
  outfile << "[" << N/100 << ", " << 2 <<"]\n";
  for(unsigned i=0;i<N-1;i++){
    // Velocity Verlet method
    // /* Comment to use the velocity Verlet algorithm

//    cout <<rEarth[1] << h <<vEarth[1] << hh << " "<<FES<<" "<<mEi<<endl;
    rEarth[2] = rEarth[1] + h*vEarth[1] + hh*0.5*FES*mEi;
  //  rSun[2] = rSun[1] + h*vSun[1] - hh*0.5*FES*mSi;

    diff = rEarth[2]-rSun[2];
    R = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    FESnew = -GMoE*diff/pow(R,3) * (1+gen_rel);

    vEarth[2] = vEarth[1] + h*0.5*(FES+FESnew)*mEi;
  //  vSun[2] = vSun[1] - h*0.5*(FES+FESnew)*mSi;

    FES = FESnew; //*/

    // Euler-Cromer method
    /*  Comment to use the Euler-Cromer algorithm
      diff = rEarth[1]-rSun[1];
      R = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
      FES = -GMoE*diff/pow(R,3) * (1+gen_rel);

      vEarth[2] = vEarth[1] + h*FES*mEi;
  //    vSun[2] = vSun[1] - h*FES*mSi;
      rEarth[2] = rEarth[1] + h*vEarth[2];
    //  rSun[2] = rSun[1] + h*vSun[1]; //*/

    // Find perihelion:
    if((rEarth[0].length() > rEarth[1].length())&&(rEarth[1].length() < rEarth[2].length())){
      cout << h*i << rEarth[1] << " = " << atan(rEarth[1][1]/rEarth[1][0])*180/M_PI*3600 << endl;
    }
    // Write to file:
    if(i%100==0){
    outfile << h*i <<";"<<rSun[1]<<" "<<vSun[1]<<" "<<rEarth[1] <<" "<< vEarth[1]<<"\n";
  }
  rEarth[0] = rEarth[1];
  vEarth[0] = vEarth[1];
  rEarth[1] = rEarth[2];
  vEarth[1] = vEarth[2];
  rSun[0] = rSun[1];
  vSun[0] = vSun[1];
  rSun[1] = rSun[2];
  vSun[1] = vSun[2];

  }
  outfile.close();

  return 0;
}
