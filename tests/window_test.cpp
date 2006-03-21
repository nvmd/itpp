/*!
 * \file
 * \brief Windowing functions test program
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger and Adam Piatyszek
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
#include <itpp/itsignal.h>

using namespace itpp;
using namespace std;

int main(void)
{

 cout << "================================" << endl;
 cout << "    Test of window functions    " << endl;
 cout << "================================" << endl;

 cout << "hamming(32) = " << hamming(32) << endl;
 cout << "hamming(128) = " << hamming(128) << endl;

 cout << "hanning(32) = " << hanning(32) << endl;
 cout << "hanning(128) = " << hanning(128) << endl;

 cout << "hann(32) = " << hann(32) << endl;
 cout << "hann(128) = " << hann(128) << endl;

 cout << "blackman(32) = " << blackman(32) << endl;
 cout << "blackman(128) = " << blackman(128) << endl;

 cout << "triang(32) = " << triang(32) << endl;
 cout << "triang(128) = " << triang(128) << endl;

 return 0;

}
