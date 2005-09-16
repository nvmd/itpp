/*!
 * \file 
 * \brief Definitions of some specific functions useful in communications
 * \author Tony Ottosson
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
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
 * -------------------------------------------------------------------------
 */

#ifndef COMMFUNC_H
#define COMMFUNC_H

#include <itpp/base/mat.h>
#include <itpp/base/vec.h>

namespace itpp {

  /*!
    \brief Generate gray code of blocklength m.
    \ingroup miscfunc

    The codes are contained as binary codewords {0,1} in the rows of the 
		returned matrix.
    See also the \c gray() function in \c math/scalfunc.h.
  */
  bmat graycode(int m);

  /*! 
    \brief Calculate the Hamming distance between \a a and \a b
    \ingroup miscfunc
  */
  int hamming_distance(const bvec &a, const bvec &b);

  /*! 
    \brief Calculate the Hamming weight of \a a
    \ingroup miscfunc
  */
  int weight(const bvec &a);

} // namespace itpp

#endif // #ifndef COMMFUNC_H
