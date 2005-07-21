/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 2005 by Johan Bergman.                                      *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Functions for Fix, Fixed, CFix and CFixed
  \author Johan Bergman
  
  $Revision$
  
  $Date$
*/

#include <itpp/fixedpoint/fix_functions.h>
#include <itpp/base/itassert.h>

namespace itpp {

  vec to_vec(const fixvec &v)
  {
    vec temp(v.length());
    for (int i=0; i<v.length(); i++) {
      temp(i) = v(i).unfix();
    }
    return temp;
  }

  cvec to_cvec(const cfixvec &v)
  {
    cvec temp(v.length());
    for (int i=0; i<v.length(); i++) {
      temp(i) = v(i).unfix();
    }
    return temp;
  }

  mat to_mat(const fixmat &m)
  {
    mat temp(m.rows(), m.cols());
    for (int i=0; i<m.rows(); i++) {
      for (int j=0; j<m.cols(); j++) {
        temp(i,j) = m(i,j).unfix();
      }
    }
    return temp;
  }

  cmat to_cmat(const cfixmat &m)
  {
    cmat temp(m.rows(), m.cols());
    for (int i=0; i<m.rows(); i++) {
      for (int j=0; j<m.cols(); j++) {
        temp(i,j) = m(i,j).unfix();
      }
    }
    return temp;
  }

  Fix abs(const Fix &x)
  {
    fixrep tmp = x.get_re();
    return Fix((tmp >= 0 ? tmp : -tmp),  // Risk for overflow!
               x.get_shift(),
               0, 0);
  }

  Fix real(const CFix &x)
  {
    return Fix(x.get_re(),
               x.get_shift(),
               0, 0);
  }

  Fix imag(const CFix &x)
  {
    return Fix(x.get_im(),
               x.get_shift(),
               0, 0);
  }

  CFix conj(const CFix &x)
  {
    return CFix(x.get_re(),
                -x.get_im(),
                x.get_shift(),
                0, 0);
  }

} //namespace itpp
