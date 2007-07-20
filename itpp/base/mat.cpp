/*!
 * \file
 * \brief Matrix Class Implementation
 * \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/mat.h>

#if defined (HAVE_BLAS)
#  include <itpp/base/blas.h>
#endif


namespace itpp {

  template<>
  bool Mat<int>::set(const char *values)
  {
    std::istringstream buffer(values);
    int rows = 0, maxrows = 10, cols = 0, nocols = 0, maxcols = 10;
    bool comma = true;
    std::streamoff offset = 0;

    free();
    alloc(maxrows, maxcols);
    zeros();

    // go through the whole buffer
    while (buffer.peek() != EOF) {
      if (++rows > maxrows) {
	maxrows <<= 1;
	set_size(maxrows, maxcols, true);
      }

      cols = 0;
      // parsing of each row starts here
      while (buffer.peek() != ';' && buffer.peek() != EOF) {
	switch (buffer.peek()) {
	  // first remove value separators
	case ',':
	  it_assert(!comma, "Mat<Num_T>::set(): Improper matrix string");
	  buffer.seekg(1, std::ios_base::cur);
	  comma = true;
	  break;
	case ' ': case '\t':
	  buffer.seekg(1, std::ios_base::cur);
	  break;

	  // then check for a possible sign
	case '+': case '-':
	  buffer.seekg(1, std::ios_base::cur);
	  offset--;
	  break;

	  // check for hexadecimal, octal or zero number
	case '0':
	  buffer.seekg(1, std::ios_base::cur);
	  offset--;
	  if (buffer.peek() != EOF && buffer.peek() != ';') {
	    switch (buffer.peek()) {
	      // hexadecimal number
	    case 'x': case 'X':
	      buffer.seekg(1, std::ios_base::cur);
	      offset--;
	      do {
		switch (buffer.peek()) {
		case '0': case '1': case '2': case '3': case '4': case '5':
		case '6': case '7': case '8': case '9': case 'a': case 'A':
		case 'b': case 'B': case 'c': case 'C': case 'd': case 'D':
		case 'e': case 'E': case 'f': case 'F':
		  buffer.seekg(1, std::ios_base::cur);
		  offset--;
		  break;
		default:
		  it_error("Mat<int>::set(): Improper data string (hex)");
		}
	      } while (buffer.peek() != EOF && buffer.peek() != ';' &&
		       buffer.peek() != ' ' && buffer.peek() != '\t' &&
		       buffer.peek() != ',');

	      buffer.clear();
	      buffer.seekg(offset, std::ios_base::cur);
	      offset = 0;
	      if (++cols > nocols) {
		// it_assert(rows == 1, "Mat<int>::set(): Improper number of columns in matrix string");
		nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	      }
	      buffer >> std::hex >> this->operator()(rows-1, cols-1);
	      comma = false;
	      break; // end of hexadeximal part

	      // octal number
	    case '1': case '2': case '3': case '4':
	    case '5': case '6': case '7':
	      buffer.seekg(1, std::ios_base::cur);
	      offset--;
	      while (buffer.peek() != EOF && buffer.peek() != ';' &&
		     buffer.peek() != ' ' && buffer.peek() != '\t' &&
		     buffer.peek() != ',') {
		switch (buffer.peek()) {
		case '0': case '1': case '2': case '3':
		case '4': case '5': case '6': case '7':
		  buffer.seekg(1, std::ios_base::cur);
		  offset--;
		  break;
		default:
		  it_error("Mat<int>::set(): Improper data string (oct)");
		}
	      }
	      buffer.clear();
	      buffer.seekg(offset, std::ios_base::cur);
	      offset = 0;
	      if (++cols > nocols) {
		// it_assert(rows == 1, "Mat<int>::set(): Improper number of columns in matrix string");
		nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	      }
	      buffer >> std::oct >> this->operator()(rows-1, cols-1);
	      comma = false;
	      break; // end of octal part

	      // zero
	    case ' ': case '\t': case ',': case ';':
	      buffer.clear();
	      buffer.seekg(offset, std::ios_base::cur);
	      offset = 0;
	      if (++cols > nocols) {
		// it_assert(rows == 1, "Mat<int>::set(): Improper number of columns in matrix string");
		nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	      }
	      buffer >> std::dec >> this->operator()(rows-1, cols-1);
	      comma = false;
	      break; // end of zero part

	    default:
	      it_error("Mat<int>::set(): Improper data string");
	    }
	  } // if (buffer.peek() != EOF && buffer.peek() != ';')
	  else { // parse zero at the end of buffer/row
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    if (++cols > nocols) {
	      // it_assert(rows == 1, "Mat<int>::set(): Improper number of columns in matrix string");
	      nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	    }
	    buffer >> std::dec >> this->operator()(rows-1, cols-1);
	    comma = false;
	  }
	  break; // end of check for hex, oct or zero number

	  // decimal number
	case '1': case '2': case '3': case '4': case '5':
	case '6': case '7': case '8': case '9':
	  buffer.seekg(1, std::ios_base::cur);
	  offset--;
	  while (buffer.peek() != EOF && buffer.peek() != ';' &&
		 buffer.peek() != ' ' && buffer.peek() != '\t' &&
		 buffer.peek() != ',') {
	    switch (buffer.peek()) {
	    case '0': case '1': case '2': case '3': case '4':
	    case '5': case '6': case '7': case '8': case '9':
	      buffer.seekg(1, std::ios_base::cur);
	      offset--;
	      break;
	    default:
	      it_error("Mat<int>::set(): Improper data string (dec)");
	    }
	  }
	  buffer.clear();
	  buffer.seekg(offset, std::ios_base::cur);
	  offset = 0;
	  if (++cols > nocols) {
	    // it_assert(rows == 1, "Mat<int>::set(): Improper number of columns in matrix string");
	    nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	  }
	  buffer >> std::dec >> this->operator()(rows-1, cols-1);
	  comma = false;
	  break; // end of decimal part

	default:
	  it_error("Mat<int>::set(): Improper data string");
	}
      } // while (buffer.peek() != ';' && buffer.peek() != EOF)

      // check if the current row has the same number of cols as the first row
      // it_assert(cols == nocols, "Mat<int>::set(): Improper number of cols in matrix string");

      if (buffer.peek() == ';') {
	// remove row separator...
	buffer.seekg(1, std::ios_base::cur);
	// ... but also check if it was at the end of buffer
	while (buffer.peek() == ' ' || buffer.peek() == '\t') {
	  buffer.seekg(1, std::ios_base::cur);
	}
	comma = false;
	// it_assert(buffer.peek() != EOF, "Mat<int>::set(): Improper data string");
      }

      it_assert(!comma || (nocols == 0), "Mat<int>::set(): Improper matrix string");

    } // while (buffer.peek() != EOF)

    set_size(rows, nocols, true);

    return true;
  }

  template<>
  bool Mat<short int>::set(const char *values) {
    std::istringstream buffer(values);
    int rows = 0, maxrows = 10, cols = 0, nocols = 0, maxcols = 10;
    bool comma = true;
    std::streamoff offset = 0;

    free();
    alloc(maxrows, maxcols);
    zeros();

    // go through the whole buffer
    while (buffer.peek() != EOF) {
      if (++rows > maxrows) {
	maxrows <<= 1;
	set_size(maxrows, maxcols, true);
      }

      cols = 0;
      // parsing of each row starts here
      while (buffer.peek() != ';' && buffer.peek() != EOF) {
	switch (buffer.peek()) {
	  // first remove value separators
	case ',':
	  it_assert(!comma, "Mat<short>::set(): Improper matrix string");
	  buffer.seekg(1, std::ios_base::cur);
	  comma = true;
	  break;
	case ' ': case '\t':
	  buffer.seekg(1, std::ios_base::cur);
	  break;

	  // then check for a possible sign
	case '+': case '-':
	  buffer.seekg(1, std::ios_base::cur);
	  offset--;
	  break;

	  // check for hexadecimal, octal or zero number
	case '0':
	  buffer.seekg(1, std::ios_base::cur);
	  offset--;
	  if (buffer.peek() != EOF && buffer.peek() != ';') {
	    switch (buffer.peek()) {
	      // hexadecimal number
	    case 'x': case 'X':
	      buffer.seekg(1, std::ios_base::cur);
	      offset--;
	      do {
		switch (buffer.peek()) {
		case '0': case '1': case '2': case '3': case '4': case '5':
		case '6': case '7': case '8': case '9': case 'a': case 'A':
		case 'b': case 'B': case 'c': case 'C': case 'd': case 'D':
		case 'e': case 'E': case 'f': case 'F':
		  buffer.seekg(1, std::ios_base::cur);
		  offset--;
		  break;
		default:
		  it_error("Mat<short int>::set(): Improper data string (hex)");
		}
	      } while (buffer.peek() != EOF && buffer.peek() != ';' &&
		       buffer.peek() != ' ' && buffer.peek() != '\t' &&
		       buffer.peek() != ',');

	      buffer.clear();
	      buffer.seekg(offset, std::ios_base::cur);
	      offset = 0;
	      if (++cols > nocols) {
		// it_assert(rows == 1, "Mat<short int>::set(): Improper number of columns in matrix string");
		nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	      }
	      buffer >> std::hex >> this->operator()(rows-1, cols-1);
	      comma = false;
	      break; // end of hexadeximal part

	      // octal number
	    case '1': case '2': case '3': case '4':
	    case '5': case '6': case '7':
	      buffer.seekg(1, std::ios_base::cur);
	      offset--;
	      while (buffer.peek() != EOF && buffer.peek() != ';' &&
		     buffer.peek() != ' ' && buffer.peek() != '\t' &&
		     buffer.peek() != ',') {
		switch (buffer.peek()) {
		case '0': case '1': case '2': case '3':
		case '4': case '5': case '6': case '7':
		  buffer.seekg(1, std::ios_base::cur);
		  offset--;
		  break;
		default:
		  it_error("Mat<short int>::set(): Improper data string (oct)");
		}
	      }
	      buffer.clear();
	      buffer.seekg(offset, std::ios_base::cur);
	      offset = 0;
	      if (++cols > nocols) {
		// it_assert(rows == 1, "Mat<short int>::set(): Improper number of columns in matrix string");
		nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	      }
	      buffer >> std::oct >> this->operator()(rows-1, cols-1);
	      comma = false;
	      break; // end of octal part

	      // zero
	    case ' ': case '\t': case ',': case ';':
	      buffer.clear();
	      buffer.seekg(offset, std::ios_base::cur);
	      offset = 0;
	      if (++cols > nocols) {
		// it_assert(rows == 1, "Mat<short int>::set(): Improper number of columns in matrix string");
		nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	      }
	      buffer >> std::dec >> this->operator()(rows-1, cols-1);
	      comma = false;
	      break; // end of zero part

	    default:
	      it_error("Mat<short int>::set(): Improper data string");
	    }
	  } // if (buffer.peek() != EOF && buffer.peek() != ';')
	  else { // parse zero at the end of buffer/row
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    if (++cols > nocols) {
	      // it_assert(rows == 1, "Mat<short int>::set(): Improper number of columns in matrix string");
	      nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	    }
	    buffer >> std::dec >> this->operator()(rows-1, cols-1);
	    comma = false;
	  }
	  break; // end of check for hex, oct or zero number

	  // decimal number
	case '1': case '2': case '3': case '4': case '5':
	case '6': case '7': case '8': case '9':
	  buffer.seekg(1, std::ios_base::cur);
	  offset--;
	  while (buffer.peek() != EOF && buffer.peek() != ';' &&
		 buffer.peek() != ' ' && buffer.peek() != '\t' &&
		 buffer.peek() != ',') {
	    switch (buffer.peek()) {
	    case '0': case '1': case '2': case '3': case '4':
	    case '5': case '6': case '7': case '8': case '9':
	      buffer.seekg(1, std::ios_base::cur);
	      offset--;
	      break;
	    default:
	      it_error("Mat<short int>::set(): Improper data string (dec)");
	    }
	  }
	  buffer.clear();
	  buffer.seekg(offset, std::ios_base::cur);
	  offset = 0;
	  if (++cols > nocols) {
	    // it_assert(rows == 1, "Mat<short int>::set(): Improper number of columns in matrix string");
	    nocols = cols;
		if (cols > maxcols) {
		  maxcols <<= 1;
		  set_size(maxrows, maxcols, true);
		}
	  }
	  buffer >> std::dec >> this->operator()(rows-1, cols-1);
	  comma = false;
	  break; // end of decimal part

	default:
	  it_error("Mat<short int>::set(): Improper data string");
	}
      } // while (buffer.peek() != ';' && buffer.peek() != EOF)

      // check if the current row has the same number of cols as the first row
      // it_assert(cols == nocols, "Mat<short int>::set(): Improper number of cols in matrix string");

      if (buffer.peek() == ';') {
	// remove row separator...
	buffer.seekg(1, std::ios_base::cur);
	// ... but also check if it was at the end of buffer
	while (buffer.peek() == ' ' || buffer.peek() == '\t') {
	  buffer.seekg(1, std::ios_base::cur);
	}
	comma = false;
	// it_assert(buffer.peek() != EOF, "Mat<short int>::set(): Improper data string");
      }
      it_assert(!comma || (nocols == 0), "Mat<short>::set(): Improper matrix string");

    } // while (buffer.peek() != EOF)

    set_size(rows, nocols, true);

    return true;
  }


  template<>
  const cmat Mat<std::complex<double> >::hermitian_transpose() const
  {
    cmat temp(no_cols, no_rows);
    for (int i=0; i<no_rows; i++)
      for (int j=0; j<no_cols; j++)
	temp(j,i) = std::conj(operator()(i,j));

    return temp;
  }


  // -------- Multiplication operator -------------

#if defined(HAVE_BLAS)

  template<>
  mat& Mat<double>::operator*=(const mat &m)
  {
    it_assert_debug(no_cols == m.no_rows,"mat::operator*=: wrong sizes");
    mat r(no_rows, m.no_cols); // unnecessary memory??
    double alpha = 1.0;
    double beta = 0.0;
    char trans = 'n';
    blas::dgemm_(&trans, &trans, &no_rows, &m.no_cols, &no_cols, &alpha, data,
		 &no_rows, m.data, &m.no_rows, &beta, r.data, &r.no_rows);
    operator=(r); // time consuming
    return *this;
  }

  template<>
  cmat& Mat<std::complex<double> >::operator*=(const cmat &m)
  {
    it_assert_debug(no_cols == m.no_rows,"cmat::operator*=: wrong sizes");
    cmat r(no_rows, m.no_cols); // unnecessary memory??
    std::complex<double> alpha = std::complex<double>(1.0);
    std::complex<double> beta = std::complex<double>(0.0);
    char trans = 'n';
    blas::zgemm_(&trans, &trans, &no_rows, &m.no_cols, &no_cols, &alpha, data,
		 &no_rows, m.data, &m.no_rows, &beta, r.data, &r.no_rows);
    operator=(r); // time consuming
    return *this;
  }

  template<>
  const mat operator*(const mat &m1, const mat &m2)
  {
    it_assert_debug(m1.no_cols == m2.no_rows,"mat::operator*: wrong sizes");
    mat r(m1.no_rows, m2.no_cols);
    double alpha = 1.0;
    double beta = 0.0;
    char trans = 'n';
    blas::dgemm_(&trans, &trans, &m1.no_rows, &m2.no_cols, &m1.no_cols, &alpha,
		 m1.data, &m1.no_rows, m2.data, &m2.no_rows, &beta, r.data,
		 &r.no_rows);
    return r;
  }

  template<>
  const cmat operator*(const cmat &m1, const cmat &m2)
  {
    it_assert_debug(m1.no_cols == m2.no_rows,"cmat::operator*: wrong sizes");
    cmat r(m1.no_rows, m2.no_cols);
    std::complex<double> alpha = std::complex<double>(1.0);
    std::complex<double> beta = std::complex<double>(0.0);
    char trans = 'n';
    blas::zgemm_(&trans, &trans, &m1.no_rows, &m2.no_cols, &m1.no_cols, &alpha,
		 m1.data, &m1.no_rows, m2.data, &m2.no_rows, &beta, r.data,
		 &r.no_rows);
    return r;
  }

  template<>
  const Vec<double> operator*(const mat &m, const Vec<double> &v)
  {
    it_assert_debug(m.no_cols == v.size(), "mat::operator*: wrong sizes");
    vec r(m.no_rows);
    double alpha = 1.0;
    double beta = 0.0;
    char trans = 'n';
    int incr = 1;
    blas::dgemv_(&trans, &m.no_rows, &m.no_cols, &alpha, m.data, &m.no_rows,
		 v._data(), &incr, &beta, r._data(), &incr);
    return r;
  }

  template<>
  const Vec<std::complex<double> > operator*(const cmat &m, const Vec<std::complex<double> > &v)
  {
    it_assert_debug(m.no_cols == v.size(), "cmat::operator*: wrong sizes");
    cvec r(m.no_rows);
    std::complex<double> alpha = std::complex<double>(1.0);
    std::complex<double> beta = std::complex<double>(0.0);
    char trans = 'n';
    int incr = 1;
    blas::zgemv_(&trans, &m.no_rows, &m.no_cols, &alpha, m.data, &m.no_rows,
		 v._data(), &incr, &beta, r._data(), &incr);
    return r;
  }

#endif // HAVE_BLAS

  // ---------------------- Instantiations --------------------------------

  //------- class instantiations --------

  template class Mat<double>;
  template class Mat<std::complex<double> >;
  template class Mat<int>;
  template class Mat<short int>;
  template class Mat<bin>;

  //-------- Addition operators ---------------

  template const mat operator+(const mat &m1, const mat &m2);
  template const cmat operator+(const cmat &m1, const cmat &m2);
  template const imat operator+(const imat &m1, const imat &m2);
  template const smat operator+(const smat &m1, const smat &m2);
  template const bmat operator+(const bmat &m1, const bmat &m2);

  template const mat operator+(const mat &m, double t);
  template const cmat operator+(const cmat &m, std::complex<double> t);
  template const imat operator+(const imat &m, int t);
  template const smat operator+(const smat &m, short t);
  template const bmat operator+(const bmat &m, bin t);

  template const mat operator+(double t, const mat &m);
  template const cmat operator+(std::complex<double> t, const cmat &m);
  template const imat operator+(int t, const imat &m);
  template const smat operator+(short t, const smat &m);
  template const bmat operator+(bin t, const bmat &m);

  //-------- Subraction operators ---------------

  template const mat operator-(const mat &m1, const mat &m2);
  template const cmat operator-(const cmat &m1, const cmat &m2);
  template const imat operator-(const imat &m1, const imat &m2);
  template const smat operator-(const smat &m1, const smat &m2);
  template const bmat operator-(const bmat &m1, const bmat &m2);

  template const mat operator-(const mat &m, double t);
  template const cmat operator-(const cmat &m, std::complex<double> t);
  template const imat operator-(const imat &m, int t);
  template const smat operator-(const smat &m, short t);
  template const bmat operator-(const bmat &m, bin t);

  template const mat operator-(double t, const mat &m);
  template const cmat operator-(std::complex<double> t, const cmat &m);
  template const imat operator-(int t, const imat &m);
  template const smat operator-(short t, const smat &m);
  template const bmat operator-(bin t, const bmat &m);

  //--------- Unary minus ---------------

  template const mat operator-(const mat &m);
  template const cmat operator-(const cmat &m);
  template const imat operator-(const imat &m);
  template const smat operator-(const smat &m);
  template const bmat operator-(const bmat &m);

  //-------- Multiplication operators ---------------

#if !defined(HAVE_BLAS)
  template const mat operator*(const mat &m1, const mat &m2);
  template const cmat operator*(const cmat &m1, const cmat &m2);
#endif
  template const imat operator*(const imat &m1, const imat &m2);
  template const smat operator*(const smat &m1, const smat &m2);
  template const bmat operator*(const bmat &m1, const bmat &m2);

#if !defined(HAVE_BLAS)
  template const vec operator*(const mat &m, const vec &v);
  template const cvec operator*(const cmat &m, const cvec &v);
#endif
  template const ivec operator*(const imat &m, const ivec &v);
  template const svec operator*(const smat &m, const svec &v);
  template const bvec operator*(const bmat &m, const bvec &v);

  template const mat operator*(const vec &v, const mat &m);
  template const cmat operator*(const cvec &v, const cmat &m);
  template const imat operator*(const ivec &v, const imat &m);
  template const smat operator*(const svec &v, const smat &m);
  template const bmat operator*(const bvec &v, const bmat &m);

  template const mat operator*(const mat &m, double t);
  template const cmat operator*(const cmat &m, std::complex<double> t);
  template const imat operator*(const imat &m, int t);
  template const smat operator*(const smat &m, short t);
  template const bmat operator*(const bmat &m, bin t);

  template const mat operator*(double t, const mat &m);
  template const cmat operator*(std::complex<double> t, const cmat &m);
  template const imat operator*(int t, const imat &m);
  template const smat operator*(short t, const smat &m);
  template const bmat operator*(bin t, const bmat &m);

  // ------------ Elementwise multiplication -----------

  template const mat elem_mult(const mat &A, const mat &B);
  template const cmat elem_mult(const cmat &A, const cmat &B);
  template const imat elem_mult(const imat &A, const imat &B);
  template const smat elem_mult(const smat &A, const smat &B);
  template const bmat elem_mult(const bmat &A, const bmat &B);

  template void elem_mult_out(const mat &A, const mat &B, mat &out);
  template void elem_mult_out(const cmat &A, const cmat &B, cmat &out);
  template void elem_mult_out(const imat &A, const imat &B, imat &out);
  template void elem_mult_out(const smat &A, const smat &B, smat &out);
  template void elem_mult_out(const bmat &A, const bmat &B, bmat &out);

  template void elem_mult_out(const mat &A, const mat &B, const mat &C, mat &out);
  template void elem_mult_out(const cmat &A, const cmat &B, const cmat &C, cmat &out);
  template void elem_mult_out(const imat &A, const imat &B, const imat &C, imat &out);
  template void elem_mult_out(const smat &A, const smat &B, const smat &C, smat &out);
  template void elem_mult_out(const bmat &A, const bmat &B, const bmat &C, bmat &out);

  template void elem_mult_out(const mat &A, const mat &B, const mat &C, const mat &D, mat &out);
  template void elem_mult_out(const cmat &A, const cmat &B, const cmat &C, const cmat &D, cmat &out);
  template void elem_mult_out(const imat &A, const imat &B, const imat &C, const imat &D, imat &out);
  template void elem_mult_out(const smat &A, const smat &B, const smat &C, const smat &D, smat &out);
  template void elem_mult_out(const bmat &A, const bmat &B, const bmat &C, const bmat &D, bmat &out);

  template void elem_mult_inplace(const mat &A, mat &B);
  template void elem_mult_inplace(const cmat &A, cmat &B);
  template void elem_mult_inplace(const imat &A, imat &B);
  template void elem_mult_inplace(const smat &A, smat &B);
  template void elem_mult_inplace(const bmat &A, bmat &B);

  template double elem_mult_sum(const mat &A, const mat &B);
  template std::complex<double> elem_mult_sum(const cmat &A, const cmat &B);
  template int elem_mult_sum(const imat &A, const imat &B);
  template short elem_mult_sum(const smat &A, const smat &B);
  template bin elem_mult_sum(const bmat &A, const bmat &B);

  // ------------ Division operator -----------

  template const mat operator/(const mat &m, double t);
  template const cmat operator/(const cmat &m, std::complex<double> t);
  template const imat operator/(const imat &m, int t);
  template const smat operator/(const smat &m, short t);
  template const bmat operator/(const bmat &m, bin t);

  // ------------ Elementwise division -----------

  template const mat elem_div(const mat &A, const mat &B);
  template const cmat elem_div(const cmat &A, const cmat &B);
  template const imat elem_div(const imat &A, const imat &B);
  template const smat elem_div(const smat &A, const smat &B);
  template const bmat elem_div(const bmat &A, const bmat &B);

  template void elem_div_out(const mat &A, const mat &B, mat &out);
  template void elem_div_out(const cmat &A, const cmat &B, cmat &out);
  template void elem_div_out(const imat &A, const imat &B, imat &out);
  template void elem_div_out(const smat &A, const smat &B, smat &out);
  template void elem_div_out(const bmat &A, const bmat &B, bmat &out);

  template double elem_div_sum(const mat &A, const mat &B);
  template std::complex<double> elem_div_sum(const cmat &A, const cmat &B);
  template int elem_div_sum(const imat &A, const imat &B);
  template short elem_div_sum(const smat &A, const smat &B);
  template bin elem_div_sum(const bmat &A, const bmat &B);

  // ------------- Concatenations -----------------

  template const mat concat_horizontal(const mat &m1, const mat &m2);
  template const cmat concat_horizontal(const cmat &m1, const cmat &m2);
  template const imat concat_horizontal(const imat &m1, const imat &m2);
  template const smat concat_horizontal(const smat &m1, const smat &m2);
  template const bmat concat_horizontal(const bmat &m1, const bmat &m2);

  template const mat concat_vertical(const mat &m1, const mat &m2);
  template const cmat concat_vertical(const cmat &m1, const cmat &m2);
  template const imat concat_vertical(const imat &m1, const imat &m2);
  template const smat concat_vertical(const smat &m1, const smat &m2);
  template const bmat concat_vertical(const bmat &m1, const bmat &m2);

  //----------- Output stream --------------

  template std::ostream &operator<<(std::ostream &os, const mat  &m);
  template std::ostream &operator<<(std::ostream &os, const cmat &m);
  template std::ostream &operator<<(std::ostream &os, const imat  &m);
  template std::ostream &operator<<(std::ostream &os, const smat  &m);
  template std::ostream &operator<<(std::ostream &os, const bmat  &m);

  //----------- Input stream -------------

  template std::istream &operator>>(std::istream &is, mat  &m);
  template std::istream &operator>>(std::istream &is, cmat &m);
  template std::istream &operator>>(std::istream &is, imat  &m);
  template std::istream &operator>>(std::istream &is, smat  &m);
  template std::istream &operator>>(std::istream &is, bmat  &m);


} // namespace itpp
