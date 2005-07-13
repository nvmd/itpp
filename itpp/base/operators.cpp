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
  \brief Implementation of operators for vectors and matricies of different types.
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include "itpp/base/operators.h"
#include "itpp/base/binary.h"

namespace itpp {

  //-----------  Scalar and a ivec -----------------

  vec operator+(const double &s, const ivec &v)
  {
    it_assert1(v.size() > 0, "operator+(): Vector of zero length");

    vec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = s + double(v(i));
    }
    return temp;
  }

  vec operator-(const double &s, const ivec &v)
  {
    it_assert1(v.size() > 0, "operator-(): Vector of zero length");

    vec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = s - double(v(i));
    }
    return temp;
  }

  vec operator*(const double &s, const ivec &v)
  {
    it_assert1(v.size() > 0, "operator*(): Vector of zero length");

    vec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = s * double(v(i));
    }
    return temp;
  }

  vec operator/(const double &s, const ivec &v)
  {
    it_assert1(v.size() > 0, "operator/(): Vector of zero length");

    vec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = s / double(v(i));
    }
    return temp;
  }

  vec operator/(const ivec &v, const double &s)
  {
    it_assert1(v.size() > 0, "operator/(): Vector of zero length");

    vec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = double(v(i))/s;
    }
    return temp;
  }

  cvec operator+(const std::complex<double> &s, const ivec &v)
  {
    it_assert1(v.size() > 0, "operator+(): Vector of zero length");

    cvec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = s + std::complex<double>(v(i));
    }
    return temp;
  }

  cvec operator-(const std::complex<double> &s, const ivec &v)
  {
    it_assert1(v.size() > 0, "operator-(): Vector of zero length");

    cvec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = s - std::complex<double>(v(i));
    }
    return temp;
  }

  cvec operator*(const std::complex<double> &s, const ivec &v)
  {
    it_assert1(v.size() > 0, "operator*(): Vector of zero length");

    cvec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = s * std::complex<double>(v(i));
    }
    return temp;
  }

  cvec operator/(const std::complex<double> &s, const ivec &v)
  {
    it_assert1(v.size() > 0, "operator/(): Vector of zero length");

    cvec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = s / std::complex<double>(v(i));
    }
    return temp;
  }

  cvec operator/(const ivec &v, const std::complex<double> &s)
  {
    it_assert1(v.size() > 0, "operator/(): Vector of zero length");

    cvec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = std::complex<double>(v(i))/s;
    }
    return temp;
  }

  //-----------  Scalar and a cvec -----------------

  cvec operator+(const double &s, const cvec &v)
  {
    it_assert1(v.size() > 0, "operator+(): Vector of zero length");

    cvec temp = v;
    for (int i=0;i<v.size();i++) {
      temp(i) += s;
    }
    return temp;
  }

  cvec operator-(const double &s, const cvec &v)
  {
    it_assert1(v.size() > 0, "operator-(): Vector of zero length");

    cvec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i) = std::complex<double>((double)s - v(i).real(), -v(i).imag());
    }
    return temp;
  }

  cvec operator*(const double &s, const cvec &v)
  {
    it_assert1(v.size() > 0, "operator*(): Vector of zero length");

    cvec temp = v;
    for (int i=0;i<v.size();i++) {
      temp(i) *= (double)s; 
    }
    return temp;
  }

  cvec operator/(const cvec &v, const double &s)
  {
    it_assert1(v.size() > 0, "operator/(): Vector of zero length");

    cvec temp = v;
    for (int i=0;i<v.size();i++) {
      temp(i) /= (double)s;
    }
    return temp;
  }

  cvec operator/(const double &s, const cvec &v)
  {
    it_assert1(v.size() > 0, "operator/(): Vector of zero length");

    cvec temp(v.length());
    for (int i=0;i<v.size();i++) {
      temp(i) = s / v(i);
    }
    return temp;
  }

  //-----------  Scalar and a cmat -----------------

  cmat operator+(const double &s, const cmat &m)
  {
    it_assert1(m.rows()> 0 && m.cols() > 0, "operator+(): Matrix of zero length");

    cmat temp = m;
    for (int i=0;i<m._datasize();i++) {
      temp._data()[i] += s;
    }
    return temp;
  }

  cmat operator-(const double &s, const cmat &m)
  {
    it_assert1(m.rows()> 0 && m.cols() > 0, "operator-(): Matrix of zero length");

    cmat temp(m.rows(), m.cols());
    for (int i=0;i<m._datasize();i++) {
      temp._data()[i] = std::complex<double>((double)s - m(i).real(), -m(i).imag());
    }
    return temp;
  }

  cmat operator*(const double &s, const cmat &m)
  {
    it_assert1(m.rows()> 0 && m.cols() > 0, "operator*(): Matrix of zero length");

    cmat temp = m;
    for (int i=0;i<m._datasize();i++) {
      temp._data()[i] *= (double)s; 
    }
    return temp;
  }

  cmat operator*(const std::complex<double> &s, const mat &m)
  {
    it_assert1(m.rows()> 0 && m.cols() > 0, "operator*(): Matrix of zero length");

    cmat temp(m.rows(), m.cols());

    for (int i=0;i<m._datasize();i++) {
      temp._data()[i] = s*m._data()[i]; 
    }
    return temp;
  }

  cmat operator/(const cmat &m, const double &s)
  {
    it_assert1(m.rows()> 0 && m.cols() > 0, "operator/(): Matrix of zero length");

    cmat temp = m;
    for (int i=0;i<m._datasize();i++) {
      temp._data()[i] /= (double)s;
    }
    return temp;
  }

  //---------------------- between matrix and scalar --------------------

  //----------- Multiplication of a scalar and a vec -----------------

  cvec operator*(const std::complex<double> &s, const vec &v)
  {
    cvec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i)=s*std::complex<double>(v(i),0.0);
    }
    return temp;
  }

  cvec operator*(const vec &v, const std::complex<double> &s)
  {
    cvec temp(v.size());
    for (int i=0;i<v.size();i++) {
      temp(i)=s*std::complex<double>(v(i),0.0);
    }
    return temp;
  }

  // ===============================================================================================

  // ---------------- Addition of vectors ---------------


  vec operator+(const bvec &a, const vec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes does not match");
    vec temp(a.size());
    for (int i=0;i<a.size();i++) {temp(i)=(double)a(i)+b(i);}
    return temp;
  }

  vec operator+(const svec &a, const vec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes does not match");
    vec temp(a.size());
    for (int i=0;i<a.size();i++) {temp(i)=(double)a(i)+b(i);}
    return temp;
  }

  vec operator+(const ivec &a, const vec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes does not match");
    vec temp(a.size());
    for (int i=0;i<a.size();i++) {temp(i)=(double)a(i)+b(i);}
    return temp;
  }

  cvec operator+(const bvec &a, const cvec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes does not match");
    cvec temp = b;
    for (int i=0;i<a.size();i++) {temp(i) +=(double)a(i);}
    return temp;
  }

  cvec operator+(const svec &a, const cvec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes does not match");
    cvec temp = b;
    for (int i=0;i<a.size();i++) {temp(i) +=(double)a(i);}
    return temp;
  }

  cvec operator+(const ivec &a, const cvec &b)
  {
    it_assert1(a.size() == b.size(), "operator+(): sizes does not match");
    cvec temp = b;
    for (int i=0;i<a.size();i++) {temp(i) +=(double)a(i);}
    return temp;
  }

  // ---------------- Multiplication of vectors ---------------

  double operator*(const bvec &a, const vec &b)
  {
    it_assert1(a.size() == b.size(), "operator*(): sizes does not match");
    double temp=0;
    for (int i=0;i<a.size();i++) {temp+=(double)a(i)*b(i);}
    return temp;
  }

  double operator*(const svec &a, const vec &b)
  {
    it_assert1(a.size() == b.size(), "operator*(): sizes does not match");
    double temp=0;
    for (int i=0;i<a.size();i++) {temp+=(double)a(i)*b(i);}
    return temp;
  }

  double operator*(const ivec &a, const vec &b)
  {
    it_assert1(a.size() == b.size(), "operator*(): sizes does not match");
    double temp=0;
    for (int i=0;i<a.size();i++) {temp+=(double)a(i)*b(i);}
    return temp;
  }

  std::complex<double> operator*(const bvec &a, const cvec &b)
  {
    it_assert1(a.size() == b.size(), "operator*(): sizes does not match");
    std::complex<double> temp=0;
    for (int i=0;i<a.size();i++) {temp+=(double)a(i)*b(i);}
    return temp;
  }

  std::complex<double> operator*(const svec &a, const cvec &b)
  {
    it_assert1(a.size() == b.size(), "operator*(): sizes does not match");
    std::complex<double> temp=0;
    for (int i=0;i<a.size();i++) {temp+=(double)a(i)*b(i);}
    return temp;
  }

  std::complex<double> operator*(const ivec &a, const cvec &b)
  {
    it_assert1(a.size() == b.size(), "operator*(): sizes does not match");
    std::complex<double> temp=0;
    for (int i=0;i<a.size();i++) {temp+=(double)a(i)*b(i);}
    return temp;
  }

  // ---------------- Addition of matricies ---------------

  mat operator+(const bmat &a, const mat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes does not match");
    mat temp(b);

    for (int i=0;i<a.rows();i++) {
      for (int j=0;j<a.cols();j++) {
	temp(i,j)+=(double)a(i,j);
      }
    }
    return temp;
  }

  mat operator+(const smat &a, const mat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes does not match");
    mat temp(b);

    for (int i=0;i<a.rows();i++) {
      for (int j=0;j<a.cols();j++) {
	temp(i,j)+=(double)a(i,j);
      }
    }
    return temp;
  }

  mat operator+(const imat &a, const mat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes does not match");
    mat temp(b);

    for (int i=0;i<a.rows();i++) {
      for (int j=0;j<a.cols();j++) {
	temp(i,j)+=(double)a(i,j);
      }
    }
    return temp;
  }

  // ---------------- Addition of cmat and matrices ---------------

  cmat operator+(const bmat &a, const cmat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes does not match");
    cmat temp(b);

    for (int i=0;i<a.rows();i++) {
      for (int j=0;j<a.cols();j++) {
	temp(i,j)+=(std::complex<double>)a(i,j);
      }
    }
    return temp;
  }

  cmat operator+(const smat &a, const cmat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes does not match");
    cmat temp(b);

    for (int i=0;i<a.rows();i++) {
      for (int j=0;j<a.cols();j++) {
	temp(i,j)+=(double)a(i,j);
      }
    }
    return temp;
  }

  cmat operator+(const imat &a, const cmat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes does not match");
    cmat temp(b);

    for (int i=0;i<a.rows();i++) {
      for (int j=0;j<a.cols();j++) {
	temp(i,j)+=(std::complex<double>)a(i,j);
      }
    }
    return temp;
  }

  cmat operator+(const mat &a, const cmat &b)
  {
    it_assert1(a.cols()==b.cols() && a.rows()==b.rows(), "operator+(): sizes does not match");
    cmat temp(b);

    for (int i=0;i<a.rows();i++) {
      for (int j=0;j<a.cols();j++) {
	temp(i,j)+=(std::complex<double>)a(i,j);
      }
    }
    return temp;
  }


} //namespace itpp
