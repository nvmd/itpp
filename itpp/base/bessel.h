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
  \brief Definitions of Bessel functions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __bessel_h
#define __bessel_h

#include "itpp/base/vec.h"

namespace itpp {

  /*! \addtogroup besselfunctions
   */

  /*! 
    \ingroup besselfunctions
    \brief Bessel function of first kind of order \a nu for \a nu integer

    The bessel function of first kind is defined as:
    \f[
    J_{\nu}(x) = \sum_{k=0}^{\infty} \frac{ (-1)^{k} }{k! \Gamma(\nu+k+1) } \left(\frac{x}{2}\right)^{\nu+2k}
    \f]
    where \f$\nu\f$ is the order and \f$ 0 < x < \infty \f$.
  */
  double besselj(int nu, double x);

  /*! 
    \ingroup besselfunctions
    \brief Bessel function of first kind of order \a nu for \a nu integer
  */
  vec besselj(int nu, const vec &x);

  /*! 
    \ingroup besselfunctions
    \brief Bessel function of first kind of order \a nu. \a nu is real.
  */
  double besselj(double nu, double x);

  /*! 
    \ingroup besselfunctions
    \brief Bessel function of first kind of order \a nu. \a nu is real.
  */
  vec besselj(double nu, const vec &x);

  /*! 
    \ingroup besselfunctions
    \brief Bessel function of second kind of order \a nu. \a nu is integer.

    The Bessel function of second kind is defined as:
    \f[
    Y_{\nu}(x) = \frac{J_{\nu}(x) \cos(\nu\pi) - J_{-\nu}(x)}{\sin(\nu\pi)} 
    \f]
    where \f$\nu\f$ is the order and \f$ 0 < x < \infty \f$.
  */
  double bessely(int nu, double x);

  /*! 
    \ingroup besselfunctions
    \brief Bessel function of second kind of order \a nu. \a nu is integer.
  */
  vec bessely(int nu, const vec &x);

  /*!   
    \ingroup besselfunctions
    \brief Bessel function of second kind of order \a nu. \a nu is real.
  */
  double bessely(double nu, double x);

  /*!   
    \ingroup besselfunctions
    \brief Bessel function of second kind of order \a nu. \a nu is real.
  */
  vec bessely(double nu, const vec &x);

  /*! 
    \ingroup besselfunctions
    \brief Modified Bessel function of first kind of order \a nu. \a nu is \a double. \a x is \a double.

    The Modified Bessel function of first kind is defined as:
    \f[
    I_{\nu}(x) = i^{-\nu} J_{\nu}(ix)
    \f]
    where \f$\nu\f$ is the order and \f$ 0 < x < \infty \f$.
  */
  double besseli(double nu, double x);

  /*! 
    \ingroup besselfunctions
    \brief Modified Bessel function of first kind of order \a nu. \a nu is \a double. \a x is \a double.
  */
  vec besseli(double nu, const vec &x);

  /*! 
    \ingroup besselfunctions
    \brief Modified Bessel function of second kind of order \a nu. \a nu is double. \a x is double.

    The Modified Bessel function of second kind is defined as:
    \f[
    K_{\nu}(x) = \frac{\pi}{2} i^{\nu+1} [J_{\nu}(ix) + i Y_{\nu}(ix)]
    \f]
    where \f$\nu\f$ is the order and \f$ 0 < x < \infty \f$.
  */
  double besselk(int nu, double x);

  /*! 
    \ingroup besselfunctions
    \brief Modified Bessel function of second kind of order \a nu. \a nu is double. \a x is double.
  */
  vec besselk(int nu, const vec &x);

} //namespace itpp

#endif // __bessel_h
