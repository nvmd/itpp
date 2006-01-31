/*!
* \file
* \brief Test program for the LLR class
* \author Erik G. Larsson
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


#include <itpp/itcomm.h>

using std::cout;
using std::endl;
using namespace itpp;

int main()
{
 LLR_calc_unit lcu1;
 LLR_calc_unit lcu2(10,7,9);  // low table resolution
 LLR_calc_unit lcu3(2,15,0);  // low table resolution and low LLR granuality 
 LLR_calc_unit lcu4(10,0,0);  // this gives logexp=logmax
 cout << lcu1 << endl;
 cout << lcu2 << endl;
 cout << lcu3 << endl;

 cout << "Testing Jacobian logarithm with four different resolutions." << endl;
 for (double x=0.0; x<10; x+=0.1) {
   cout << "JacLog(" << x << ") = " 
	 << lcu1.to_double(lcu1.logexp(lcu1.to_qllr(x))) << " ; " 
	 << lcu2.to_double(lcu2.logexp(lcu2.to_qllr(x))) << " ; " 
	 << lcu3.to_double(lcu3.logexp(lcu3.to_qllr(x))) << " ; " 
	 << lcu4.to_double(lcu4.logexp(lcu4.to_qllr(x)))  
	 << endl;
 }

 cout << "-------------------" << endl;
 cout << "Some special cases:" << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(100.0),lcu1.to_qllr(0.75))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(100.0),lcu1.to_qllr(-0.75))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-100.0),lcu1.to_qllr(0.75))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-100.0),lcu1.to_qllr(-0.75))) << endl;

 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(0.0),lcu1.to_qllr(0.75))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(0.0),lcu1.to_qllr(-0.75))) << endl;

 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(3.75),lcu1.to_qllr(-1.25))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-1.25),lcu1.to_qllr(3.75))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(-3.75),lcu1.to_qllr(1.25))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(1.25),lcu1.to_qllr(-3.75))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(3.75),lcu1.to_qllr(1.25))) << endl;
 cout << lcu1.to_double(lcu1.Boxplus(lcu1.to_qllr(1.25),lcu1.to_qllr(3.75))) << endl;
 
}

