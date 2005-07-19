//
// source.cpp
//
// This program is supposed to test the source classes
//
// $Id$
//

#include <iostream>
#include <iomanip>

#include "itpp/base/matfunc.h"
#include "itpp/base/random.h"
#include "itpp/base/source.h"
#include "itpp/base/stat.h"

using std::cout;
using std::endl;
using namespace itpp;

#define N 100000

using namespace std;

void show(const char *name, double sm, double sv)
{
  cout << setw(20) << name << "  "
       << setw(12) << sm << "  "
       << setw(12) << sv << endl;
}

void show(const char *name, complex<double> sm, double sv)
{
  cout << setw(20) << name << "  "
       << setw(12) << sm << "  "
       << setw(12) << sv << endl;
}

#define REALRUN(name,s)           \
for (i=0; i<N; i++)               \
    real_result(i) = s();              \
show(name,mean(real_result),variance(real_result));

#define COMPLEXRUN(name,s)        \
for (i=0; i<N; i++)               \
    complex_result(i) = s();              \
show(name,mean(complex_result),variance(complex_result));

int main()
{
  int i;
  vec real_result;
  cvec complex_result;

  Uniform_RNG        s0;
  Exponential_RNG    s1;
  Normal_RNG         s2;
  Weibull_RNG        s3;
  Rayleigh_RNG       s4;
  I_Uniform_RNG      s5(3, 8);
  AR1_Normal_RNG     s6(2.0, 1.0, 0.95);
  Complex_Normal_RNG s7;

  Sine_Source      s10(20.0/N);
  Square_Source    s11(20.0/N);
  Triangle_Source  s12(20.0/N);
  Sawtooth_Source  s13(20.0/N);
  Impulse_Source   s14(20.0/N);
  Pattern_Source   s15(vec("1 3"));

  RNG_reset(12345);

  cout << setw(20) << "Source" << "  "
       << setw(12) << "sim mean" << "  "
       << setw(12) << "sim var" << endl;
  for (int i=0; i<20+4*14; i++)
    cout << '=';
  cout << endl;

  real_result.set_size(N,false);
  complex_result.set_size(N,false);

  REALRUN("Uniform",        s0);
  REALRUN("Exponential",    s1);
  REALRUN("Normal",         s2);
  REALRUN("Weibull",        s3);
  REALRUN("Rayleigh",       s4);
  REALRUN("I_Uniform",      s5);
  REALRUN("AR1_Normal",     s6);
  COMPLEXRUN("Complex_Normal", s7);

  REALRUN("Sine",        s10);
  REALRUN("Square",      s11);
  REALRUN("Triangle",    s12);
  REALRUN("Sawtooth",    s13);
  REALRUN("Impulse",     s14);
  REALRUN("Pattern",     s15);

  return 0;
}
