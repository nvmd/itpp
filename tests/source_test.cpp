/*!
 * \file 
 * \brief Deterministic sources test program
 * \author Tobias Ringstrom, Tony Ottosson and Adam Piatyszek
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
 */

#include <iomanip>
#include <itpp/itbase.h>

using namespace itpp;
using namespace std;

#define LOOP_SIZE 100000
#define THRESHOLD 1e-13

#define REALRUN(name,s)					\
  for (i = 0; i < LOOP_SIZE; i++)			\
    real_result(i) = s();				\
  show(name, mean(real_result), variance(real_result));

#define COMPLEXRUN(name,s)					\
  for (i = 0; i < LOOP_SIZE; i++)				\
    complex_result(i) = s();					\
  show(name, mean(complex_result), variance(complex_result));


void show(const char *name, double sm, double sv)
{
  cout << setw(18) << name << "  "
       << setw(20) << round_to_zero(sm, THRESHOLD) << "  "
       << setw(20) << round_to_zero(sv, THRESHOLD) << endl;
}


void show(const char *name, complex<double> sm, double sv)
{
  cout << setw(18) << name << "  "
       << setw(20) << round_to_zero(sm, THRESHOLD) << "  "
       << setw(20) << round_to_zero(sv, THRESHOLD) << endl;
}


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

  Sine_Source      s10(20.0/LOOP_SIZE);
  Square_Source    s11(20.0/LOOP_SIZE);
  Triangle_Source  s12(20.0/LOOP_SIZE);
  Sawtooth_Source  s13(20.0/LOOP_SIZE);
  Impulse_Source   s14(20.0/LOOP_SIZE);
  Pattern_Source   s15(vec("1 3"));

  RNG_reset(12345);

  cout << setw(18) << "Source" << "  "
       << setw(20) << "sim mean" << "  "
       << setw(20) << "sim var" << endl
       << "============================================================================" << endl;

  real_result.set_size(LOOP_SIZE, false);
  complex_result.set_size(LOOP_SIZE, false);

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
