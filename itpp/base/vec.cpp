/*!
 * \file
 * \brief Templated Vector Class Implementation
 * \author Tony Ottosson and Tobias Ringstrom
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

#include <itpp/base/vec.h>

#if defined (HAVE_CBLAS) || defined(HAVE_MKL)
#include <itpp/base/cblas.h>
#endif

namespace itpp {

  template<>
  bool cvec::set(const char *values)
  {
    std::istringstream buffer(values);
    int pos=0, maxpos=10;

    alloc(maxpos);
    
    while (buffer.peek() != EOF) {

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
				while (buffer.peek() == ' ') { buffer.get(); }
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

    while (buffer.peek() != EOF) {
      if (buffer.peek() == ',') {
				buffer.get();
      } else {
				pos++;
				if (pos > maxpos) {
					maxpos *= 2;
					set_size(maxpos, true);
				}
				buffer >> intemp;
				data[pos-1] = intemp;
				while (buffer.peek() == ' ') { buffer.get(); }
			}
    }
    set_size(pos, true);

    return true;
  }

#if defined(HAVE_CBLAS) || defined(HAVE_MKL)
  template<>
  const double dot(const vec &v1, const vec &v2)
  {
    it_assert1(v1.datasize == v2.datasize, "vec::dot: wrong sizes");
    double r=0.0;

    r= cblas_ddot(v1.datasize, v1.data, 1, v2.data, 1);

    return r;
  }
 
  template<>
  const std::complex<double> dot(const cvec &v1, const cvec &v2)
  {
    it_assert1(v1.datasize == v2.datasize, "cvec::dot: wrong sizes");
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

  template const vec operator+(const vec &v1, const vec &v2);
  template const cvec operator+(const cvec &v1, const cvec &v2);
  template const ivec operator+(const ivec &v1, const ivec &v2);
  template const svec operator+(const svec &v1, const svec &v2);
  template const bvec operator+(const bvec &v1, const bvec &v2);

  template const vec operator+(const vec &v1, double t);
  template const cvec operator+(const cvec &v1, std::complex<double> t);
  template const ivec operator+(const ivec &v1, int t);
  template const svec operator+(const svec &v1, short t);
  template const bvec operator+(const bvec &v1, bin t);

  template const vec operator+(double t, const vec &v1);
  template const cvec operator+(std::complex<double> t, const cvec &v1);
  template const ivec operator+(int t, const ivec &v1);
  template const svec operator+(short t, const svec &v1);
  template const bvec operator+(bin t, const bvec &v1);

  //------------- Subraction operator ----------

  template const vec operator-(const vec &v1, const vec &v2);
  template const cvec operator-(const cvec &v1, const cvec &v2);
  template const ivec operator-(const ivec &v1, const ivec &v2);
  template const svec operator-(const svec &v1, const svec &v2);
  template const bvec operator-(const bvec &v1, const bvec &v2);

  template const vec operator-(const vec &v, double t);
  template const cvec operator-(const cvec &v, std::complex<double> t);
  template const ivec operator-(const ivec &v, int t);
  template const svec operator-(const svec &v, short t);
  template const bvec operator-(const bvec &v, bin t);

  template const vec operator-(double t, const vec &v);
  template const cvec operator-(std::complex<double> t, const cvec &v);
  template const ivec operator-(int t, const ivec &v);
  template const svec operator-(short t, const svec &v);
  template const bvec operator-(bin t, const bvec &v);

  //---------- Unary minus -------------

  template const vec operator-(const vec &v);
  template const cvec operator-(const cvec &v);
  template const ivec operator-(const ivec &v);
  template const svec operator-(const svec &v);
  template const bvec operator-(const bvec &v);

  //------------- Multiplication operator ----------

  template const double dot(const vec &v1, const vec &v2);
  template const std::complex<double> dot(const cvec &v1, const cvec &v2);
  template const int dot(const ivec &v1, const ivec &v2);
  template const short dot(const svec &v1, const svec &v2);
  template const bin dot(const bvec &v1, const bvec &v2);

  template const int operator*(const ivec &v1, const ivec &v2);
  template const short operator*(const svec &v1, const svec &v2);
  template const bin operator*(const bvec &v1, const bvec &v2);

  template const mat outer_product(const vec &v1, const vec &v2);
  template const cmat outer_product(const cvec &v1, const cvec &v2);
  template const imat outer_product(const ivec &v1, const ivec &v2);
  template const smat outer_product(const svec &v1, const svec &v2);
  template const bmat outer_product(const bvec &v1, const bvec &v2);

  template const vec operator*(const vec &v, double t);
  template const cvec operator*(const cvec &v, std::complex<double> t);
  template const ivec operator*(const ivec &v, int t);
  template const svec operator*(const svec &v, short t);
  template const bvec operator*(const bvec &v, bin t);

  template const vec operator*(double t, const vec &v);
  template const cvec operator*(std::complex<double> t, const cvec &v);
  template const ivec operator*(int t, const ivec &v);
  template const svec operator*(short t, const svec &v);
  template const bvec operator*(bin t, const bvec &v);

  //------------- Elementwise Multiplication operator (two vectors) ----------

  template const vec elem_mult(const vec &v1, const vec &v2);
  template const cvec elem_mult(const cvec &v1, const cvec &v2);
  template const ivec elem_mult(const ivec &v1, const ivec &v2);
  template const svec elem_mult(const svec &v1, const svec &v2);
  template const bvec elem_mult(const bvec &v1, const bvec &v2);

  //------------- Elementwise Multiplication operator (three vectors) ----------

  template const vec elem_mult(const vec &v1, const vec &v2, const vec &v3);
  template const cvec elem_mult(const cvec &v1, const cvec &v2, const cvec &v3);
  template const ivec elem_mult(const ivec &v1, const ivec &v2, const ivec &v3);
  template const svec elem_mult(const svec &v1, const svec &v2, const svec &v3);
  template const bvec elem_mult(const bvec &v1, const bvec &v2, const bvec &v3);

  //------------- Elementwise Multiplication operator (four vectors) ----------

  template const vec elem_mult(const vec &v1, const vec &v2, const vec &v3, const vec &v4);
  template const cvec elem_mult(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4);
  template const ivec elem_mult(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4);
  template const svec elem_mult(const svec &v1, const svec &v2, const svec &v3, const svec &v4);
  template const bvec elem_mult(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4);

  //------------- Division operator ----------

  template const vec operator/(const vec &v, double t);
  template const cvec operator/(const cvec &v, std::complex<double> t);
  template const ivec operator/(const ivec &v, int t);
  template const svec operator/(const svec &v, short t);
  template const bvec operator/(const bvec &v, bin t);

  template const vec operator/(const double t, const vec &v);
  template const cvec operator/(const std::complex<double> t, const cvec &v);
  template const ivec operator/(const int t, const ivec &v);
  template const svec operator/(const short t, const svec &v);
  template const bvec operator/(const bin t, const bvec &v);

  //------------- Elementwise Division operator ----------

  template const vec elem_div(const vec &v1, const vec &v2);
  template const cvec elem_div(const cvec &v1, const cvec &v2);
  template const ivec elem_div(const ivec &v1, const ivec &v2);
  template const svec elem_div(const svec &v1, const svec &v2);
  template const bvec elem_div(const bvec &v1, const bvec &v2);

  template const vec elem_div(const double t, const vec &v);
  template const cvec elem_div(const std::complex<double> t, const cvec &v);
  template const ivec elem_div(const int t, const ivec &v);
  template const svec elem_div(const short t, const svec &v);
  template const bvec elem_div(const bin t, const bvec &v);

  //--------------------- concat operator -----------------

  template const vec concat(const vec &v, const double a);
  template const cvec concat(const cvec &v, const std::complex<double> a);
  template const ivec concat(const ivec &v, const int a);
  template const svec concat(const svec &v, const short a);
  template const bvec concat(const bvec &v, const bin a);

  template const vec concat(const double a, const vec &v);
  template const cvec concat(const std::complex<double> a, const cvec &v);
  template const ivec concat(const int a, const ivec &v);
  template const svec concat(const short a, const svec &v);
  template const bvec concat(const bin a, const bvec &v);

  template const vec concat(const vec &v1, const vec &v2);
  template const cvec concat(const cvec &v1, const cvec &v2);
  template const ivec concat(const ivec &v1, const ivec &v2);
  template const svec concat(const svec &v1, const svec &v2);
  template const bvec concat(const bvec &v1, const bvec &v2);

  template const vec concat(const vec &v1, const vec &v2, const vec &v3);
  template const cvec concat(const cvec &v1, const cvec &v2, const cvec &v3);
  template const ivec concat(const ivec &v1, const ivec &v2, const ivec &v3);
  template const svec concat(const svec &v1, const svec &v2, const svec &v3);
  template const bvec concat(const bvec &v1, const bvec &v2, const bvec &v3);

  template const vec concat(const vec &v1, const vec &v2, const vec &v3, const vec &v4);
  template const cvec concat(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4);
  template const ivec concat(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4);
  template const svec concat(const svec &v1, const svec &v2, const svec &v3, const svec &v4);
  template const bvec concat(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4);

  template const vec concat(const vec &v1, const vec &v2, const vec &v3, const vec &v4, const vec &v5);
  template const cvec concat(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4, const cvec &v5);
  template const ivec concat(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4, const ivec &v5);
  template const svec concat(const svec &v1, const svec &v2, const svec &v3, const svec &v4, const svec &v5);
  template const bvec concat(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4, const bvec &v5);

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

} // namespace itpp
