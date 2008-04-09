/*!
 * \file
 * \brief IT++ compatibility types and functions
 * \author Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2008  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#include <itpp/base/itcompat.h>
#include <limits>

//! \cond

#ifndef HAVE_EXPM1
#ifndef M_LN2
#define M_LN2 0.69314718055994530941723212146 // ln(2)
#endif
// Implementation taken from GSL (http://www.gnu.org/software/gsl/)
double expm1(double x)
{
  /* FIXME: this should be improved */
  if (std::fabs(x) < M_LN2) {
    /* Compute the taylor series S = x + (1/2!) x^2 + (1/3!) x^3 + ... */
    double i = 1.0;
    double sum = x;
    double term = x / 1.0;
    do {
      i++;
      term *= x / i;
      sum += term;
    }
    while (std::fabs(term) > (std::fabs(sum)
                              * std::numeric_limits<double>::epsilon()));
    return sum;
  }
  else {
    return std::exp(x) - 1.0;
  }
}
#endif // HAVE_EXPM1

//! \endcond

