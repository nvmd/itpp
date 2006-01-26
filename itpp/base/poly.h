/*!
 * \file
 * \brief Polynomial functions
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

#ifndef POLY_H
#define POLY_H

#include <itpp/base/vec.h>


namespace itpp {

  /*!
    \addtogroup poly
  */



  /*!
    \brief Create a polynomial of the given roots
    \ingroup poly
    
    Create a polynomial \c p with roots \c r
  */
  //!@{
  void poly(const vec &r, vec &p);
  inline vec poly(const vec &r) { vec temp; poly(r, temp); return temp; }
  void poly(const cvec &r, cvec &p);
  inline cvec poly(const cvec &r) { cvec temp; poly(r, temp); return temp; }
  //!@}
  

  /*!
    \brief Calculate the roots of the polynomial
    \ingroup poly
    
    Calculate the roots \c r of the polynomial \c p
  */
  //!@{
  void roots(const vec &p, cvec &r);
  inline cvec roots(const vec &p) { cvec temp; roots(p, temp); return temp; }
  void roots(const cvec &p, cvec &r);
  inline cvec roots(const cvec &p) { cvec temp; roots(p, temp); return temp; }
  //!@}
  

  /*!
    \brief Evaluate polynomial
    \ingroup poly
    
    Evaluate the polynomial \c p (of length \f$N+1\f$ at the points \c x
    The output is given by
    \f[
    p_0 x^N + p_1 x^{N-1} + \ldots + p_{N-1} x + p_N
    \f]
  */  
  //!@{
  vec polyval(const vec &p, const vec &x);
  cvec polyval(const vec &p, const cvec &x);
  cvec polyval(const cvec &p, const vec &x);
  cvec polyval(const cvec &p, const cvec &x);
  //!@}

} // namespace itpp

#endif // #ifndef POLY_H
