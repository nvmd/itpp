/*!
 * \file 
 * \brief Implementation of some specific functions useful in communications
 * \author Tony Ottosson
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

#include <itpp/comm/commfunc.h>
#include <itpp/base/converters.h>
#include <itpp/base/specmat.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/binary.h>


namespace itpp { 

  bmat graycode(int m)
  {
    if (m == 1) {
      smat temp = "0;1";
      return to_bmat(temp);
    } else {
      bvec temp(1<<(m-1));
      bmat bb	= graycode(m-1);
      bmat out(1<<m, m);
      out.zeros();
      out.set_col(0, concat(zeros_b(1<<(m-1)), ones_b(1<<(m-1))) );
      for (int i=0; i<m-1; i++) {
				temp = bb.get_col(i);
				out.set_col(i+1, concat(temp, reverse(temp)) );
      }
      return out;
    }
  }

  int hamming_distance(const bvec &a, const bvec &b)
  {
    int i, n=0;

    it_assert1(a.size() == b.size(), "hamming_distance()");
    for (i=0; i<a.size(); i++)
      if (a(i) != b(i))
				n++;
    
    return n;
  }

  int weight(const bvec &a)
  {
    int i, n=0;

    for (i=0; i<a.size(); i++)
      if (a(i)==bin(1))
				n++;
    
    return n;
  }

} // namespace itpp
