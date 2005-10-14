/*!
 * \file
 * \brief Help functions to make functions with vec and mat as arguments
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

#ifndef HELP_FUNCTIONS_H
#define HELP_FUNCTIONS_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>

namespace itpp {

  /*! 
    \brief Help function to call for a function: vec function(vec)
    \ingroup miscfunc
  */
  inline vec vec_function(double (*f)(double), const vec &x)
    {
      vec out(x.length());
      for (int i=0;i<x.length();i++) 
	out(i) = f(x(i));
      return out;
    }

  /*! 
    \brief Help function to call for a function: cvec function(cvec)
    \ingroup miscfunc
  */
  inline cvec cvec_function(std::complex<double> (*f)(const std::complex<double> &), const cvec &x)
    {
      cvec out(x.length());
      for (int i=0;i<x.length();i++) 
	out(i) = f(x(i));
      return out;
    }

  /*! 
    \brief Help function to call for a function: mat function(mat)
    \ingroup miscfunc
  */
  inline mat mat_function(double (*f)(double), const mat &x)
    {
      mat out(x.rows(),x.cols());
      for (int i=0;i<out.rows();i++) {
	for (int j=0;j<out.cols();j++) {
	  out(i,j) = f(x(i,j));
	}
      }
      return out;
    }
  /*! 
    \brief Help function to call for a function: cmat function(cmat)
    \ingroup miscfunc
  */

  inline cmat cmat_function(std::complex<double> (*f)(const std::complex<double> &), const cmat &x)
    {
      cmat out(x.rows(),x.cols());
      for (int i=0;i<out.rows();i++) {
	for (int j=0;j<out.cols();j++) {
	  out(i,j) = f(x(i,j));
	}
      }
      return out;
    }

  /*! 
    \brief Help function to call for a function: vec function(double,vec)
    \ingroup miscfunc
  */
  inline vec double_vec_function(double (*f)(double,double), const double x, const vec &y)
    {
      vec out(y.length());
      for (int i=0;i<y.length();i++) 
	out(i) = f(x,y(i));
      return out;
    }

  /*! 
    \brief Help function to call for a function: mat function(double,mat)
    \ingroup miscfunc
  */
  inline mat double_mat_function(double (*f)(double,double), const double x, const mat &y)
    {
      mat out(y.rows(),y.cols());
      for (int i=0;i<out.rows();i++) {
	for (int j=0;j<out.cols();j++) {
	  out(i,j) = f(x, y(i,j));
	}
      }
      return out;
    }

  /*! 
    \brief Help function to call for a function: vec function(vec,double)
    \ingroup miscfunc
  */
  inline vec vec_double_function(double (*f)(double,double), const vec &x, const double y)
    {
      vec out(x.length());
      for (int i=0;i<x.length();i++) 
	out(i) = f(x(i),y);
      return out;
    }

  /*! 
    \brief Help function to call for a function: mat function(mat,double)
    \ingroup miscfunc
  */
  inline mat mat_double_function(double (*f)(double,double), const mat &x, const double y)
    {
      mat out(x.rows(),x.cols());
      for (int i=0;i<out.rows();i++) {
	for (int j=0;j<out.cols();j++) {
	  out(i,j) = f(x(i,j), y);
	}
      }
      return out;
    }

} // namespace itpp

#endif // #ifndef HELP_FUNCTIONS_H

