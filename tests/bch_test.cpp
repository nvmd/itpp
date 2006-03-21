/*!
 * \file 
 * \brief BCH encoder/decoder class test program
 * \author Pal Frenger and Adam Piatyszek
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

#include <itpp/itbase.h>
#include <itpp/itcomm.h>

using namespace itpp;
using namespace std;


bvec set_errors(const bvec &input, const ivec errpos) 
{
  bvec output = input;
  for (int i = 0; i < errpos.length(); i++)
    output(errpos(i)) ^= 1;
  return output ;
}


int main() 
{
  BCH bch(31, 21, 2, "3 5 5 1");

  cout << "===================================" << endl;
  cout << "    Test of BCH encoder/decoder    " << endl;
  cout << "===================================" << endl;

  cout << "A two error case (should be corrected)" << endl;
  bvec input = randb(21);
  bvec encoded = bch.encode(input);
  bvec err = set_errors(encoded, (ivec) "1 2"); // error positions
  bvec decoded = bch.decode(err);

  cout << "input =   " << input << endl ;
  cout << "encoded = " << encoded << endl ;
  cout << "err =     " << err <<  endl ;
  cout << "decoded = " << decoded << endl ;
 
  input = randb(21);
  encoded = bch.encode(input);
  err = set_errors(encoded, (ivec) "1 2 27"); // error positions
  decoded = bch.decode(err);

  cout << "A three error case (will cause decoding errors)" << endl;
  cout << "input =   " << input << endl;
  cout << "encoded = " << encoded << endl;
  cout << "err =     " << err <<  endl;
  cout << "decoded = " << decoded << endl;
 
  return 0;
}
