/*!
* \file 
* \brief Timer classes test program
* \author Thomas Eriksson, Tony Ottosson, Tobias Ringstrom and Adam Piatyszek
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

using namespace std;
using namespace itpp;

int main()
{
   CPU_Timer t1;
   Real_Timer t2;

   tic();
   t1.start();
   t2.start();

   while (t2.get_time() <= 4.0)
     ;

   t1.stop();
   t2.stop();
   double t3 = toc();

   if (t2.get_time() >= 4.0)
     cout << "Real_Timer is OK" << endl;

   // Check against 1 ms margin
   if (t1.get_time() - t2.get_time() <= 1e-3)
     cout << "CPU_Timer is OK" << endl;

   if (t3 >= t2.get_time())
     cout << "tic() and toc() is OK" << endl;

   return 0;
}
