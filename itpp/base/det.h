/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2005 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Determinant calculation of square matrices
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __det_h
#define __det_h

#include "itpp/base/mat.h"

namespace itpp {


  /*!
    \addtogroup determinant
  */

  /*!
    \brief Determinant of real square matrix.
    \ingroup determinant

    Calculate determinant of the real matrix \f$\mathbf{X}\f$

    Uses LU-factorisation.
    \f[
    \det(\mathbf{X}) = \det(\mathbf{P}^T \mathbf{L}) \det(\mathbf{U}) = \det(\mathbf{P}^T) \prod(\mathrm{diag}(\mathbf{U}))
    \f]
    and the determinant of the permuation matrix is \f$ \pm 1\f$ dependening on the number of row permuations
  */
  double det(const mat &X);


  /*!
    \brief Determinant of complex square matrix.
    \ingroup determinant

    Calculate determinant of the complex matrix \f$\mathbf{X}\f$

    Uses LU-factorisation.
    \f[
    \det(\mathbf{X}) = \det(\mathbf{P}^T \mathbf{L}) \det(\mathbf{U}) = \det(\mathbf{P}^T) \prod(\mathrm{diag}(\mathbf{U}))
    \f]
    and the determinant of the permuation matrix is \f$ \pm 1\f$ dependening on the number of row permuations
  */
  std::complex<double> det(const cmat &X);


}//namespace itpp

#endif // __det_h

