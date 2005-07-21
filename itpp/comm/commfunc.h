/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2001 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Definitions of some specific functions useful in communications
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef _commfunc_h
#define _commfunc_h

#include <itpp/base/mat.h>
#include <itpp/base/vec.h>

namespace itpp {

  /*!
    \brief Generate gray code of blocklength m.
    \ingroup miscfunc

    The codes are contained as binary codewords {0,1} in the rows of the returned matrix.
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

} //namespace itpp

#endif // _commfunc_h
