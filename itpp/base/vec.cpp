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
  \brief Templated Vector Class Implementation
  \author Tony Ottosson and Tobias Ringström

  $Revision$

  $Date$ 
*/

#include "itpp/base/vec.h"

#if defined (HAVE_CBLAS) || defined(HAVE_MKL)
#include "itpp/base/cblas.h"
#endif

namespace itpp {

  template<>
  bool cvec::set(const char *values)
  {
    std::istringstream buffer(values);
    int pos=0, maxpos=10;

    alloc(maxpos);
    
    while (buffer.peek()!=EOF) {

      switch (buffer.peek()) {

      case ':':
        it_error("set: expressions with ':' are not valid for cvec");
	break;

      case ',':
	buffer.get();
	break;

      default:
	pos++;
	if (pos > maxpos) {
	  maxpos *= 2;
	  set_size(maxpos, true);
	}
	buffer >> data[pos-1];
	while (buffer.peek()==' ') { buffer.get(); }
	break;
      }

    }
    set_size(pos, true);

    return true;
  }

  template<>
  bool bvec::set(const char *values)
  {
    std::istringstream buffer(values);
    int pos=0, maxpos=10;
    short intemp;

    alloc(maxpos);

    while (buffer.peek()!=EOF) {
      if (buffer.peek()==',') {
	buffer.get();
      } else {
	pos++;
	if (pos > maxpos) {
	  maxpos=maxpos*2;
	  set_size(maxpos, true);
	}
	buffer >> intemp;
	data[pos-1]=intemp;
 	while (buffer.peek()==' ') { buffer.get(); }
     }
    }
    set_size(pos, true);

    return true;
  }

#if defined(HAVE_CBLAS) || defined(HAVE_MKL)
  template<>
  double dot(const vec &v1, const vec &v2)
  {
    it_assert1(v1.datasize==v2.datasize, "vec::dot: wrong sizes");
    double r=0.0;

    r= cblas_ddot(v1.datasize, v1.data, 1, v2.data, 1);

    return r;
  }
 
  template<>
  std::complex<double> dot(const cvec &v1, const cvec &v2)
  {
    it_assert1(v1.datasize==v2.datasize, "cvec::dot: wrong sizes");
    std::complex<double> r=0.0;

    cblas_zdotu_sub(v1.datasize, v1.data, 1, v2.data, 1, &r);

    return r;
  }

#endif // HAVE_CBLAS or HAVE_MKL

  template<> 
  bvec cvec::operator==(const std::complex<double>) const
  { 
    it_error("operator==: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec cvec::operator!=(const std::complex<double>) const
  { 
    it_error("operator!=: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec cvec::operator<=(const std::complex<double>) const
  { 
    it_error("operator<=: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec cvec::operator>(const std::complex<double>) const
  { 
    it_error("operator>: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec cvec::operator<(const std::complex<double>) const
  { 
    it_error("operator<: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec cvec::operator>=(const std::complex<double>) const
  { 
    it_error("operator>=: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<>
  Mat<std::complex<double> > cvec::hermitian_transpose() const
  {
    Mat<std::complex<double> > temp(1, datasize);
    for (int i=0; i<datasize; i++)
      temp(i) = std::conj(data[i]);

    return temp;
  }


  //---------------------------------------------------------------------
  // Instantiations
  //---------------------------------------------------------------------

  //--------- class instantiations -------------

  template class Vec<double>;
  template class Vec<int>;
  template class Vec<short int>;
  template class Vec<std::complex<double> >;
  template class Vec<bin>;

  //------------- Addition operator ----------

  template vec operator+(const vec &v1, const vec &v2);
  template cvec operator+(const cvec &v1, const cvec &v2);
  template ivec operator+(const ivec &v1, const ivec &v2);
  template svec operator+(const svec &v1, const svec &v2);
  template bvec operator+(const bvec &v1, const bvec &v2);

  template vec operator+(const vec &v1, double t);
  template cvec operator+(const cvec &v1, std::complex<double> t);
  template ivec operator+(const ivec &v1, int t);
  template svec operator+(const svec &v1, short t);
  template bvec operator+(const bvec &v1, bin t);

  template vec operator+(double t, const vec &v1);
  template cvec operator+(std::complex<double> t, const cvec &v1);
  template ivec operator+(int t, const ivec &v1);
  template svec operator+(short t, const svec &v1);
  template bvec operator+(bin t, const bvec &v1);

  //------------- Subraction operator ----------

  template vec operator-(const vec &v1, const vec &v2);
  template cvec operator-(const cvec &v1, const cvec &v2);
  template ivec operator-(const ivec &v1, const ivec &v2);
  template svec operator-(const svec &v1, const svec &v2);
  template bvec operator-(const bvec &v1, const bvec &v2);

  template vec operator-(const vec &v, double t);
  template cvec operator-(const cvec &v, std::complex<double> t);
  template ivec operator-(const ivec &v, int t);
  template svec operator-(const svec &v, short t);
  template bvec operator-(const bvec &v, bin t);

  template vec operator-(double t, const vec &v);
  template cvec operator-(std::complex<double> t, const cvec &v);
  template ivec operator-(int t, const ivec &v);
  template svec operator-(short t, const svec &v);
  template bvec operator-(bin t, const bvec &v);

  //---------- Unary minus -------------

  template vec operator-(const vec &v);
  template cvec operator-(const cvec &v);
  template ivec operator-(const ivec &v);
  template svec operator-(const svec &v);
  template bvec operator-(const bvec &v);

  //------------- Multiplication operator ----------

  template double dot(const vec &v1, const vec &v2);
  template std::complex<double> dot(const cvec &v1, const cvec &v2);
  template int dot(const ivec &v1, const ivec &v2);
  template short dot(const svec &v1, const svec &v2);
  template bin dot(const bvec &v1, const bvec &v2);

  template int operator*(const ivec &v1, const ivec &v2);
  template short operator*(const svec &v1, const svec &v2);
  template bin operator*(const bvec &v1, const bvec &v2);

  template mat outer_product(const vec &v1, const vec &v2);
  template cmat outer_product(const cvec &v1, const cvec &v2);
  template imat outer_product(const ivec &v1, const ivec &v2);
  template smat outer_product(const svec &v1, const svec &v2);
  template bmat outer_product(const bvec &v1, const bvec &v2);

  template vec operator*(const vec &v, double t);
  template cvec operator*(const cvec &v, std::complex<double> t);
  template ivec operator*(const ivec &v, int t);
  template svec operator*(const svec &v, short t);
  template bvec operator*(const bvec &v, bin t);

  template vec operator*(double t, const vec &v);
  template cvec operator*(std::complex<double> t, const cvec &v);
  template ivec operator*(int t, const ivec &v);
  template svec operator*(short t, const svec &v);
  template bvec operator*(bin t, const bvec &v);

  //------------- Elementwise Multiplication operator (two vectors) ----------

  template vec elem_mult(const vec &v1, const vec &v2);
  template cvec elem_mult(const cvec &v1, const cvec &v2);
  template ivec elem_mult(const ivec &v1, const ivec &v2);
  template svec elem_mult(const svec &v1, const svec &v2);
  template bvec elem_mult(const bvec &v1, const bvec &v2);

  //------------- Elementwise Multiplication operator (three vectors) ----------

  template vec elem_mult(const vec &v1, const vec &v2, const vec &v3);
  template cvec elem_mult(const cvec &v1, const cvec &v2, const cvec &v3);
  template ivec elem_mult(const ivec &v1, const ivec &v2, const ivec &v3);
  template svec elem_mult(const svec &v1, const svec &v2, const svec &v3);
  template bvec elem_mult(const bvec &v1, const bvec &v2, const bvec &v3);

  //------------- Elementwise Multiplication operator (four vectors) ----------

  template vec elem_mult(const vec &v1, const vec &v2, const vec &v3, const vec &v4);
  template cvec elem_mult(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4);
  template ivec elem_mult(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4);
  template svec elem_mult(const svec &v1, const svec &v2, const svec &v3, const svec &v4);
  template bvec elem_mult(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4);

  //------------- Division operator ----------

  template vec operator/(const vec &v, double t);
  template cvec operator/(const cvec &v, std::complex<double> t);
  template ivec operator/(const ivec &v, int t);
  template svec operator/(const svec &v, short t);
  template bvec operator/(const bvec &v, bin t);

  template vec operator/(const double t, const vec &v);
  template cvec operator/(const std::complex<double> t, const cvec &v);
  template ivec operator/(const int t, const ivec &v);
  template svec operator/(const short t, const svec &v);
  template bvec operator/(const bin t, const bvec &v);

  //------------- Elementwise Division operator ----------

  template vec elem_div(const vec &v1, const vec &v2);
  template cvec elem_div(const cvec &v1, const cvec &v2);
  template ivec elem_div(const ivec &v1, const ivec &v2);
  template svec elem_div(const svec &v1, const svec &v2);
  template bvec elem_div(const bvec &v1, const bvec &v2);

  template vec elem_div(const double t, const vec &v);
  template cvec elem_div(const std::complex<double> t, const cvec &v);
  template ivec elem_div(const int t, const ivec &v);
  template svec elem_div(const short t, const svec &v);
  template bvec elem_div(const bin t, const bvec &v);

  //--------------------- concat operator -----------------

  template vec concat(const vec &v, const double a);
  template cvec concat(const cvec &v, const std::complex<double> a);
  template ivec concat(const ivec &v, const int a);
  template svec concat(const svec &v, const short a);
  template bvec concat(const bvec &v, const bin a);

  template vec concat(const double a, const vec &v);
  template cvec concat(const std::complex<double> a, const cvec &v);
  template ivec concat(const int a, const ivec &v);
  template svec concat(const short a, const svec &v);
  template bvec concat(const bin a, const bvec &v);

  template vec concat(const vec &v1, const vec &v2);
  template cvec concat(const cvec &v1, const cvec &v2);
  template ivec concat(const ivec &v1, const ivec &v2);
  template svec concat(const svec &v1, const svec &v2);
  template bvec concat(const bvec &v1, const bvec &v2);

  template vec concat(const vec &v1, const vec &v2, const vec &v3);
  template cvec concat(const cvec &v1, const cvec &v2, const cvec &v3);
  template ivec concat(const ivec &v1, const ivec &v2, const ivec &v3);
  template svec concat(const svec &v1, const svec &v2, const svec &v3);
  template bvec concat(const bvec &v1, const bvec &v2, const bvec &v3);

  template vec concat(const vec &v1, const vec &v2, const vec &v3, const vec &v4);
  template cvec concat(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4);
  template ivec concat(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4);
  template svec concat(const svec &v1, const svec &v2, const svec &v3, const svec &v4);
  template bvec concat(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4);

  template vec concat(const vec &v1, const vec &v2, const vec &v3, const vec &v4, const vec &v5);
  template cvec concat(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4, const cvec &v5);
  template ivec concat(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4, const ivec &v5);
  template svec concat(const svec &v1, const svec &v2, const svec &v3, const svec &v4, const svec &v5);
  template bvec concat(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4, const bvec &v5);

  // -------------- output stream --------------------

  template std::ostream &operator<<(std::ostream& os, const vec &vect);
  template std::ostream &operator<<(std::ostream& os, const cvec &vect);
  template std::ostream &operator<<(std::ostream& os, const svec &vect);
  template std::ostream &operator<<(std::ostream& os, const ivec &vect);
  template std::ostream &operator<<(std::ostream& os, const bvec &vect);

  // -------------- input stream --------------------

  template std::istream &operator>>(std::istream& is, vec &vect);
  template std::istream &operator>>(std::istream& is, cvec &vect);
  template std::istream &operator>>(std::istream& is, svec &vect);
  template std::istream &operator>>(std::istream& is, ivec &vect);
  template std::istream &operator>>(std::istream& is, bvec &vect);

} //namespace itpp
