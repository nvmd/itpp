/*!
 * \file
 * \brief Random number generator test program
 * \author Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#include <itpp/itbase.h>
#include <iomanip>

using namespace itpp;
using namespace std;

int main()
{
  cout.setf(ios::fixed);
  cout.precision(2);

  RNG_reset(4357U);
  Bernoulli_RNG b_rng;
  cout << "Bernoulli_RNG:\n" << b_rng() << endl
       << b_rng(10) << endl
       << b_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Uniform_RNG u_rng;
  cout << "Uniform_RNG:\n" << u_rng() << endl
       << u_rng(10) << endl
       << u_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  I_Uniform_RNG iu_rng(0, 9);
  cout << "I_Uniform_RNG [0..9]:\n" << iu_rng() << endl
       << iu_rng(10) << endl
       << iu_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Normal_RNG n_rng;
  cout << "Normal_RNG:\n" << n_rng() << endl
       << n_rng(10) << endl
       << n_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Complex_Normal_RNG cn_rng;
  cout << "Complex_Normal_RNG:\n" << cn_rng() << endl
       << cn_rng(10) << endl
       << cn_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Exponential_RNG e_rng;
  cout << "Exponential_RNG:\n" << e_rng() << endl
       << e_rng(10) << endl
       << e_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Laplace_RNG lap_rng;
  cout << "Laplace_RNG:\n" << lap_rng() << endl
       << lap_rng(10) << endl
       << lap_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  AR1_Normal_RNG ar1_rng;
  cout << "AR1_Normal_RNG:\n" << ar1_rng() << endl
       << ar1_rng(10) << endl
       << ar1_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Weibull_RNG w_rng;
  cout << "Weibull_RNG:\n" << w_rng() << endl
       << w_rng(10) << endl
       << w_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Rayleigh_RNG ray_rng;
  cout << "Rayleigh_RNG:\n" << ray_rng() << endl
       << ray_rng(10) << endl
       << ray_rng(3, 5) << endl << endl;

  RNG_reset(4357U);
  Rice_RNG ric_rng;
  cout << "Rice_RNG:\n" << ric_rng() << endl
       << ric_rng(10) << endl
       << ric_rng(3, 5) << endl << endl;

  return 0;
}
