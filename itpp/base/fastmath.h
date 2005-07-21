/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2002 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
 * Tobias Ringström, and Jonas Samuelsson.                                   *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*! 
  \file 
  \brief Definitions of special operations on vectors and matricies optimized for speed
  \author Tony Ottosson and Tobias Ringstrom

  $Revision$

  $Date$
*/

#ifndef __fastmath_h
#define __fastmath_h

#include <itpp/base/binary.h>    // inclusion of this include made it possible to compile in M$VC; WHY??
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

} //namespace itpp

#endif // __fastmath_h















