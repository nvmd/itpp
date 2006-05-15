/*!
 * \file
 * \brief Elementary mathematical functions
 * \author Tony Ottosson and Adam Piatyszek
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

#include <itpp/base/elmatfunc.h>

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#ifndef HAVE_TGAMMA
// "True" gamma function
double tgamma(double x)
{
  double s = (2.50662827510730 + 190.9551718930764 / (x + 1)
	      - 216.8366818437280 / (x + 2) + 60.19441764023333
	      / (x + 3) - 3.08751323928546 / (x + 4) + 0.00302963870525
	      / (x + 5) - 0.00001352385959072596 / (x + 6)) / x;
  if (s < 0)
    return -exp((x + 0.5) * log(x + 5.5) - x - 5.5 + log(-s));
  else 
    return exp((x + 0.5) * log(x + 5.5) - x - 5.5 + log(s));
}
#endif

#if !defined(HAVE_LGAMMA) || (HAVE_DECL_SIGNGAM != 1)
// The sign of the Gamma function is returned in the external integer
// signgam declared in <math.h>. It is 1 when the Gamma function is positive
// or zero, -1 when it is negative. However, MinGW definition of lgamma()
// function does not use the global signgam variable.
int signgam;
// Logarithm of an absolute value of gamma function 
double lgamma(double x)
{
  double gam = tgamma(x);
  signgam = (gam < 0) ? -1 : 1;
  return log(fabs(gam));
}
#endif

#ifndef HAVE_CBRT
// Cubic root
double cbrt(double x) { return std::pow(x, 1.0/3.0); }
#endif


namespace itpp { 

  vec sqr(const cvec &data)
  {
    vec	temp(data.length());
    for (int i = 0; i < data.length(); i++)
      temp(i) = sqr(data(i));
    return temp;
  }

  mat sqr(const cmat &data)
  {
    mat	temp(data.rows(), data.cols());
    for (int i = 0; i < temp.rows(); i++) {
      for (int j = 0; j < temp.cols(); j++) {
	temp(i, j) = sqr(data(i, j));
      }
    }
    return temp;
  }

  vec abs(const cvec &data)
  {
    vec	temp(data.length());

    for (int i=0;i<data.length();i++)
      temp[i]=std::abs(data[i]);

    return temp;
  }

  mat abs(const cmat &data)
  {
    mat	temp(data.rows(),data.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	temp(i,j)=std::abs(data(i,j));
      }
    }

    return temp;
  }

  vec real(const cvec &data)
  {
    vec	temp(data.length());

    for (int i=0;i<data.length();i++)
      temp[i]=data[i].real();

    return temp;
  }

  mat real(const cmat &data)
  {
    mat	temp(data.rows(),data.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	temp(i,j)=data(i,j).real();
      }
    }

    return temp;
  }

  vec imag(const cvec &data)
  {
    vec	temp(data.length());

    for (int i=0;i<data.length();i++)
      temp[i]=data[i].imag();
    return temp;
  }

  mat imag(const cmat &data)
  {
    mat	temp(data.rows(),data.cols());
  
    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	temp(i,j)=data(i,j).imag();
      }
    }

    return temp;
  }

  vec arg(const cvec &data)
  {
    vec	temp(data.length());

    for (int i=0;i<data.length();i++)
      temp[i]=std::arg(data[i]);
	
    return temp;
  }

  mat arg(const cmat &data)
  {
    mat	temp(data.rows(),data.cols());
  
    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	temp(i,j)=std::arg(data(i,j));
      }
    }

    return temp;
  }

} // namespace itpp
