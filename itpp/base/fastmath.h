/*!
 * \file
 * \brief Definitions of special operations on vectors and matricies optimized 
 * for speed
 * \author Tony Ottosson and Tobias Ringstrom
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

#ifndef FASTMATH_H
#define FASTMATH_H

#include <itpp/base/binary.h> // inclusion of this include made it possible 
                              // to compile in M$VC; WHY??
#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/scalfunc.h>

namespace itpp {

  /*! 
    \relates Mat
    \brief Calculates m=m-v*v'*m
  */
  void sub_v_vT_m(mat &m, const vec &v);

  /*! 
    \relates Mat
    \brief Calculates m=m-m*v*v'
  */
  void sub_m_v_vT(mat &m, const vec &v);

} // namespace itpp

#endif // #ifndef FASTMATH_H
