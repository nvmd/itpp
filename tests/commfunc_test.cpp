/*!
 * \file 
 * \brief Commfunc test program
 * \author Erik G. Larsson
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

#include <itpp/itcomm.h>

using namespace std;
using namespace itpp;

int main()
{
  cout << graycode(5) << endl;

  bvec b1=randb(10);
  bvec b2=randb(10);
  cout << hamming_distance(b1,b2) << endl;

  cout << weight(b1) << endl;

  cout << waterfilling("1 2 3",0.001) << endl;
  cout << waterfilling("1 2 3",0.01) << endl;
  cout << waterfilling("1 2 3",0.1) << endl;
  cout << waterfilling("1 2 3",0.5) << endl;
  cout << waterfilling("1 2 3",1) << endl;
  cout << waterfilling("1 2 3",5) << endl;
  cout << waterfilling("1 2 3",100) << endl;
  return 0;
}
