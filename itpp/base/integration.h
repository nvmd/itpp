/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2004 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Numerical Integration
  \author Tony Ottosson

  $Revision$

  $Date$ 
*/

#ifndef __integration_h
#define __integration_h

#include <itpp/itconfig.h>
#include <itpp/base/vec.h>
#include <limits>




namespace itpp {

  /*! 
    \addtogroup integration
    \brief Numerical integration routines
  */
  
  //@{

  /*!
    1-dimensional numerical Simpson quadrature integration

    Calculate the 1-dimensional integral
    \f[
    \int_a^b f(x) dx
    \f]

    Uses an adaptive Simpson quadrature method. See [Gander] for more details. The integrand is
    specified as a function \code double f(double) \endcode.

    Example:
    \code
    #include "itpp/itbase.h"

    double f(const double x)
    {
      return x*log(x);
    }

    int main()
    {
      double res = quad( f, 1.5, 3.5);
      cout << "res = " << res << endl;

      return 0;
    }
    \endcode

    References:
    
    [Gander] Gander, W. and W. Gautschi, "Adaptive Quadrature - Revisited", BIT, Vol. 40, 2000, pp. 84-101.
    This document is also available at http:// www.inf.ethz.ch/personal/gander.
  */
  double quad(double (*f)(double), const double a, const double b, const double tol=std::numeric_limits<double>::epsilon());

  /*!
    1-dimensional numerical adaptive Lobatto quadrature integration
    
    Calculate the 1-dimensional integral
    \f[
    \int_a^b f(x) dx
    \f]

    Uses an adaptive Lobatto quadrature method. See [Gander] for more details. The integrand is
    specified as a function \code double f(double) \endcode.

    Example:
    \code
    #include "itpp/itbase.h"

    double f(const double x)
    {
      return x*log(x);
    }

    int main()
    {
      double res = quadl( f, 1.5, 3.5);
      cout << "res = " << res << endl;

      return 0;
    }
    \endcode

    References:
    
    [Gander] Gander, W. and W. Gautschi, "Adaptive Quadrature - Revisited", BIT, Vol. 40, 2000, pp. 84-101.
    This document is also available at http:// www.inf.ethz.ch/personal/gander.
  */
  double quadl(double (*f)(double), const double a, const double b, const double tol=std::numeric_limits<double>::epsilon());

  //@}
  
} // namespace itpp


#endif
