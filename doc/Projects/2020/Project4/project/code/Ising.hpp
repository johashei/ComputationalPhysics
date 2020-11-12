#ifndef Ising_HPP
#define Ising_HPP
#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <random>


using namespace std;

class Ising_2d{
private:
  double m_T,E,M;
  unsigned m_L;
  int **S;
  double w[17];
  mt19937_64 rng;
  uniform_real_distribution<double> uniform = std::uniform_real_distribution<double>(0.0,1.0);

  int wraparound(int i, int size){return (i%size + size) % size;}
  void AlignSpins();
  void RandomizeSpins();
public:
  Ising_2d(unsigned L, unsigned seed);
  ~Ising_2d();
  void Initialize(double T);
  void Metropolis(); // One Monte Carlo cycle
  double Get_Energy(){return E;}
  double Get_Magnetization(){return M;}
};

#endif
