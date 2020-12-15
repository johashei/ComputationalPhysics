#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <vector>
#include <string>
#include <unordered_map>

using namespace std;

///*
class myClass{
public:
  double X;
  void f(){cout << X << endl;}
  //void (myClass::*pf)();
  vector<void (myClass::*)()> funcs;
  void f(double x)
  {
    X = x;
    //pf = &myClass::f;
    f();
    //(this->*(pf))();
    funcs.push_back(&myClass::f);
    (this->*(funcs[0]))();
  }
};
//*/

double X;
void f();
void f(double);
vector<void (*)()> funcs;


int main(){
  myClass C;

  // Challenge : put this inside a member function
  /*  C.X = 4;
    C.pf = &myClass::f;
    C.f();
    (C.*(C.pf))();
    C.funcs.push_back(C.pf);
    (C.*(C.funcs[0]))();*/

  C.f(4);

  return 0;
}

void f(double x){
  X=x;
  funcs.push_back(&f);
}

void f(){
  cout << X << endl;
}
