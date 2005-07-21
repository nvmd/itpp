/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2003 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Implementation of elementary functions on vectors and matrices.
  \author Tony Ottosson
 
  1.5

  2003/05/22 08:55:19
*/
 
#include <itpp/base/elmatfunc.h>

using std::real;
using std::imag;
using std::arg;
using std::conj;

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

} //namespace itpp
