#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <numeric>
#include "vec3.h"

using namespace std;

int main(){

  // Set initial conditions:
  vec3 r0Sun(0,0,0);
  vec3 v0Sun(0,0,0);
  double mSun = 1;

  vec3 r0Earth(9.445675412475125E-01,3.223034192685008E-01, -1.992059498358081E-05);
  vec3 v0Earth(-5.831624954084838E-03,1.622307605140163E-02,-1.112702817686854E-06);
  v0Earth = v0Earth*365.25; // Values given in AU/day
  //vec3 r0Earth(1,0,0);
  //vec3 v0Earth(0,2*M_PI,0);
  cout << r0Earth.dot(v0Earth);

  double mEarth = 6e24/2e30;


  // Set simulation parameters
  unsigned N = 1e4;
  double T = 100; // sim time in years
  double h = T/(N-1);
  double hh = h*h;

  vec3 rEarth[N];
  vec3 rSun[N];
  vec3 vEarth[N];
  vec3 vSun[N];
  vec3 rCoM = (mSun*r0Sun + mEarth*r0Earth)/(mSun+mEarth);
  vec3 vCoM = (mSun*v0Sun + mEarth*v0Earth)/(mSun+mEarth);
  cout << rCoM << " " << vCoM << endl;

  rEarth[0] = r0Earth - rCoM;
  rSun[0] = r0Sun - rCoM;
  vEarth[0] = v0Earth - vCoM;
  vSun[0] = v0Sun - vCoM;

  // Velocity Verlet
  double GMoE = 4*M_PI*M_PI * mEarth; // ≈ G*M_sun*M_earth
  double mEi = 1/mEarth;
  double mSi = 1/mSun;

  vec3 diff = rEarth[0]-rSun[0];
  double R = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
  vec3 FES = -GMoE*diff/pow(R,3);
  vec3 FESnew;

  for(unsigned i=0;i<N-1;i++){
    rEarth[i+1] = rEarth[i] + h*vEarth[i] + hh*0.5*FES*mEi;
    rSun[i+1] = rSun[i] + h*vSun[i] - hh*0.5*FES*mSi;

    diff = rEarth[i+1]-rSun[i+1];
    R = sqrt(diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]);
    FESnew = -GMoE*diff/pow(R,3);

    vEarth[i+1] = vEarth[i] + h*0.5*(FES+FESnew)*mEi;
    vSun[i+1] = vSun[i] - h*0.5*(FES+FESnew)*mSi;

    FES = FESnew;
    //cout << rEarth[0] << " " << vEarth[0] << endl;
    //cin.get();
  }

  // Write to file:
  ofstream outfile;
  outfile.open("earth-sun.dat");
  outfile << "Earth(x y z vx vy vz) Sun(x y z vx vy vz)" << endl;
  outfile << N <<"\n";
  for(unsigned i=0;i<N;i++){
      outfile << rEarth[i][0] <<" "<< rEarth[i][1]<<" "<<rEarth[i][2]<<" "<<vEarth[i][0]<<" "<<vEarth[i][1]<<" "<<vEarth[i][2]<<" "<<rSun[i][0]<<" "<<rSun[i][1]<<" "<<rSun[i][2]<<" "<<vSun[i][0]<<" "<<vSun[i][1]<<" "<<vSun[i][2]<<"\n";
  }
  outfile.close();



  return 0;
}
