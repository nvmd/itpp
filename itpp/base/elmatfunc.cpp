/*!
 * \file
 * \brief Implementation of elementary functions on vectors and matrices
 * \author Tony Ottosson
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
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
 * -------------------------------------------------------------------------
 */

#include <itpp/base/elmatfunc.h>

namespace itpp { 

  vec sqr(const cvec &x)
  {
    vec temp(x.length());
    for (int i=0; i<temp.length(); i++)
      temp(i) = sqr(x(i));

    return temp;
  }

  mat sqr(const cmat &x)
  {
    mat temp(x.rows(), x.cols());
    for (int i=0; i<temp.rows(); i++) {
      for (int j=0; j<temp.cols(); j++) {
	temp(i,j) = sqr(x(i,j));
      }
    }

    return temp;
  }


  ivec abs(const ivec &data)
  {
    ivec temp(data.length());

    for (int i=0;i<data.length();i++)
      temp[i]=std::abs(data[i]);

    return temp;
  }

  imat abs(const imat &data)
  {
    imat temp(data.rows(),data.cols());

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
	temp(i,j)=std::abs(data(i,j));
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

  cvec conj(const cvec &data)
  {
    cvec	temp(data);

    for (int i=0;i<data.length();i++)
		temp(i)=std::conj(temp[i]);

    return temp;
  }

  cmat conj(const cmat &data)
  {
    cmat	temp(data);

    for (int i=0;i<temp.rows();i++) {
      for (int j=0;j<temp.cols();j++) {
		  temp(i,j)=std::conj(data(i,j));
      }
    }

    return temp;
  }

  bool all(const Vec<bin> &testvec)
  {
    for (int i=0; i<testvec.length(); i++) 
      if (!testvec(i)) return false;
    return true;
  }
  
  bool any(const Vec<bin> &testvec)
  {
    for (int i=0; i<testvec.length(); i++) 
      if (testvec(i)) return true;
    return false;
  } 

} // namespace itpp
