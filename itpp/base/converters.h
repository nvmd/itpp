/*!
 * \file
 * \brief Definitions of converters between different vector and matrix types
 * \author Tony Ottosson, Tobias Ringstrom, Pal Frenger and Adam Piatyszek
 *
 * -------------------------------------------------------------------------
 *
 * Copyright (C) 1995-2010  (see AUTHORS file for a list of contributors)
 *
 * This file is part of IT++ - a C++ library of mathematical, signal
 * processing, speech processing, and communications classes and functions.
 *
 * IT++ is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any
 * later version.
 *
 * IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along
 * with IT++.  If not, see <http://www.gnu.org/licenses/>.
 *
 * -------------------------------------------------------------------------
 */

#ifndef CONVERTERS_H
#define CONVERTERS_H

#include <itpp/base/help_functions.h>
#include <itpp/base/math/misc.h>


namespace itpp
{

//! \addtogroup convertfunc
//@{

// ----------------------------------------------------------------------
// Converters for vectors
// ----------------------------------------------------------------------

/*!
  \relatesalso Vec
  \brief Converts a Vec<T> to bvec
*/
template <class T>
bvec to_bvec(const Vec<T> &v)
{
  bvec temp(v.length());
  for (int i = 0; i < v.length(); ++i) {
    temp(i) = static_cast<bin>(v(i));
  }
  return temp;
}

/*!
  \relatesalso Vec
  \brief Converts a Vec<T> to svec
*/
template <class T>
svec to_svec(const Vec<T> &v)
{
  svec temp(v.length());
  for (int i = 0; i < v.length(); ++i) {
    temp(i) = static_cast<short>(v(i));
  }
  return temp;
}

/*!
  \relatesalso Vec
  \brief Converts a Vec<T> to ivec
*/
template <class T>
ivec to_ivec(const Vec<T> &v)
{
  ivec temp(v.length());
  for (int i = 0; i < v.length(); ++i) {
    temp(i) = static_cast<int>(v(i));
  }
  return temp;
}

/*!
  \relatesalso Vec
  \brief Converts a Vec<T> to vec
*/
template <class T>
vec to_vec(const Vec<T> &v)
{
  vec temp(v.length());
  for (int i = 0; i < v.length(); ++i) {
    temp(i) = static_cast<double>(v(i));
  }
  return temp;
}

/*!
  \relatesalso Vec
  \brief Converts a Vec<T> to cvec
*/
template <class T>
cvec to_cvec(const Vec<T> &v)
{
  cvec temp(v.length());
  for (int i = 0; i < v.length(); ++i) {
    temp(i) = std::complex<double>(static_cast<double>(v(i)), 0.0);
  }
  return temp;
}

//! \cond
template<> inline
cvec to_cvec(const cvec& v)
{
  return v;
}
//! \endcond

/*!
  \relatesalso Vec
  \brief Converts real and imaginary Vec<T> to cvec
*/
template <class T>
cvec to_cvec(const Vec<T> &real, const Vec<T> &imag)
{
  it_assert(real.length() == imag.length(),
            "to_cvec(): real and imaginary parts must have the same length");
  cvec temp(real.length());
  for (int i = 0; i < real.length(); ++i) {
    temp(i) = std::complex<double>(static_cast<double>(real(i)),
                                   static_cast<double>(imag(i)));
  }
  return temp;
}

/*!
  \relatesalso Vec
  \brief Converts an int to ivec
*/
ivec to_ivec(int s);

/*!
  \relatesalso Vec
  \brief Converts an double to vec
*/
vec to_vec(double s);

/*!
  \relatesalso Vec
  \brief Converts real and imaginary double to cvec
*/
cvec to_cvec(double real, double imag);

// ----------------------------------------------------------------------
// Converters for matrices
// ----------------------------------------------------------------------

/*!
  \relatesalso Mat
  \brief Converts a Mat<T> to bmat
*/
template <class T>
bmat to_bmat(const Mat<T> &m)
{
  bmat temp(m.rows(), m.cols());
  for (int i = 0; i < temp.rows(); ++i) {
    for (int j = 0; j < temp.cols(); ++j) {
      temp(i, j) = static_cast<bin>(m(i, j));
    }
  }
  return temp;
}

/*!
  \relatesalso Mat
  \brief Converts a Mat<T> to smat
*/
template <class T>
smat to_smat(const Mat<T> &m)
{
  smat temp(m.rows(), m.cols());
  for (int i = 0; i < temp.rows(); ++i) {
    for (int j = 0; j < temp.cols(); ++j) {
      temp(i, j) = static_cast<short>(m(i, j));
    }
  }
  return temp;
}

/*!
  \relatesalso Mat
  \brief Converts a Mat<T> to imat
*/
template <class T>
imat to_imat(const Mat<T> &m)
{
  imat temp(m.rows(), m.cols());
  for (int i = 0; i < temp.rows(); ++i) {
    for (int j = 0; j < temp.cols(); ++j) {
      temp(i, j) = static_cast<int>(m(i, j));
    }
  }
  return temp;
}

/*!
  \relatesalso Mat
  \brief Converts a Mat<T> to mat
*/
template <class T>
mat to_mat(const Mat<T> &m)
{
  mat temp(m.rows(), m.cols());
  for (int i = 0; i < temp.rows(); ++i) {
    for (int j = 0; j < temp.cols(); ++j) {
      temp(i, j) = static_cast<double>(m(i, j));
    }
  }
  return temp;
}

/*!
  \relatesalso Mat
  \brief Converts a Mat<T> to cmat
*/
template <class T>
cmat to_cmat(const Mat<T> &m)
{
  cmat temp(m.rows(), m.cols());
  for (int i = 0; i < temp.rows(); ++i) {
    for (int j = 0; j < temp.cols(); ++j) {
      temp(i, j) = std::complex<double>(static_cast<double>(m(i, j)), 0.0);
    }
  }
  return temp;
}

//! \cond
template<> inline
cmat to_cmat(const cmat& m)
{
  return m;
}
//! \endcond

/*!
  \relatesalso Mat
  \brief Converts real and imaginary Mat<T> to cmat
*/
template <class T>
cmat to_cmat(const Mat<T> &real, const Mat<T> &imag)
{
  it_assert_debug((real.rows() == imag.rows())
                  && (real.cols() == imag.cols()),
                  "to_cmat(): real and imag part sizes does not match");
  cmat temp(real.rows(), real.cols());
  for (int i = 0; i < temp.rows(); ++i) {
    for (int j = 0; j < temp.cols(); ++j) {
      temp(i, j) = std::complex<double>(static_cast<double>(real(i, j)),
                                        static_cast<double>(imag(i, j)));
    }
  }
  return temp;
}


/*!
  \brief Convert a decimal int \a index to bvec using \a length bits in the representation
*/
bvec dec2bin(int length, int index);

/*!
  \brief Convert a decimal int \a index to bvec. Value returned in \a v.
*/
void dec2bin(int index, bvec &v);

/*!
  \brief Convert a decimal int \a index to bvec with the first bit as MSB if \a msb_first == true
*/
bvec dec2bin(int index, bool msb_first = true);

/*!
  \brief Convert a bvec to decimal int with the first bit as MSB if \a msb_first == true
*/
int bin2dec(const bvec &inbvec, bool msb_first = true);

/*!
  \brief Convert ivec of octal form to bvec

  Converts from ivec containing {0,1,2,...,7} to bvec containing {0,1}.
  Removes zeros to the left if keepzeros = 0 (default).
  Example: oct2bin("3 5 5 1") returns {1 1 1 0 1 1 0 1 0 0 1}.
*/
bvec oct2bin(const ivec &octalindex, short keepzeros = 0);

/*!
  \brief Convert bvec to octal ivec

  Converts from  bvec containing {0,1} to ivec containing {0,1,2,...,7}.
  Adds zeros to the left if inbits.length() is not a factor of 3.
  Example: bin2oct("1 1 1 0 1 1 0 1 0 0 1") returns {3 5 5 1}.
*/
ivec bin2oct(const bvec &inbits);

//! Convert bvec to polar binary representation as ivec
ivec bin2pol(const bvec &inbvec);

//! Convert binary polar ivec to bvec
bvec pol2bin(const ivec &inpol);

//! Convert radians to degrees
inline double rad_to_deg(double x) { return (180.0 / itpp::pi * x); }
//! Convert degrees to radians
inline double deg_to_rad(double x) { return (itpp::pi / 180.0 * x); }

//! Round to nearest integer, return result in double
double round(double x);
//! Round to nearest integer
vec round(const vec &x);
//! Round to nearest integer
mat round(const mat &x);
//! Round to nearest integer
int round_i(double x);
//! Round to nearest integer and return ivec
ivec round_i(const vec &x);
//! Round to nearest integer and return imat
imat round_i(const mat &x);

//! Round to nearest upper integer
inline vec ceil(const vec &x) { return apply_function<double>(std::ceil, x); }
//! Round to nearest upper integer
inline mat ceil(const mat &x) { return apply_function<double>(std::ceil, x); }
//! The nearest larger integer
inline int ceil_i(double x) { return static_cast<int>(std::ceil(x)); }
//! Round to nearest upper integer
ivec ceil_i(const vec &x);
//! Round to nearest upper integer
imat ceil_i(const mat &x);

//! Round to nearest lower integer
inline vec floor(const vec &x) { return apply_function<double>(std::floor, x); }
//! Round to nearest lower integer
inline mat floor(const mat &x) { return apply_function<double>(std::floor, x); }
//! The nearest smaller integer
inline int floor_i(double x) { return static_cast<int>(std::floor(x)); }
//! Round to nearest lower integer
ivec floor_i(const vec &x);
//! Round to nearest lower integer
imat floor_i(const mat &x);


//! Round \a x to zero if \a abs(x) is smaller than \a threshold
inline double round_to_zero(double x, double threshold = 1e-14)
{
  return ((std::fabs(x) < threshold) ? 0.0 : x);
}

//! Round each part of \a x smaller than \a threshold to zero
inline std::complex<double> round_to_zero(const std::complex<double>& x,
    double threshold = 1e-14)
{
  return std::complex<double>(round_to_zero(x.real(), threshold),
                              round_to_zero(x.imag(), threshold));
}

//! Round each element to zero if element < threshold
inline vec round_to_zero(const vec &x, double threshold = 1e-14)
{
  return apply_function<double>(round_to_zero, x, threshold);
}

//! Round each element to zero if element < threshold
inline mat round_to_zero(const mat &x, double threshold = 1e-14)
{
  return apply_function<double>(round_to_zero, x, threshold);
}

//! Round each element to zero if element < threshold
cvec round_to_zero(const cvec &x, double threshold = 1e-14);

//! Round each element to zero if element < threshold
cmat round_to_zero(const cmat &x, double threshold = 1e-14);

//! Remove trailing digits, found after the decimal point, for numbers greater than threshold
inline double round_to_infty(const double in, const double threshold = 1e9)
{
  return (std::fabs(in)>threshold)?itpp::round(in):in;
}

//! Remove trailing digits, found after the decimal point, for complex numbers whose real and imaginary parts are greater than threshold
inline std::complex<double> round_to_infty(const std::complex<double> &in, const double threshold = 1e9)
{
  return std::complex<double>(round_to_infty(in.real(), threshold),
                              round_to_infty(in.imag(), threshold));
}

//! Remove trailing digits, found after the decimal point, for vectors greater than threshold
inline vec round_to_infty(const vec &in, const double threshold = 1e9)
{
  return apply_function<double>(round_to_infty, in, threshold);
}

//! Remove trailing digits, found after the decimal point, for matrices greater than threshold
inline mat round_to_infty(const mat &in, const double threshold = 1e9)
{
  return apply_function<double>(round_to_infty, in, threshold);
}

//! Remove trailing digits, found after the decimal point, for complex vectors greater than threshold
cvec round_to_infty(const cvec &in, const double threshold = 1e9);

//! Remove trailing digits, found after the decimal point, for complex matrices greater than threshold
cmat round_to_infty(const cmat &in, const double threshold = 1e9);

//! Convert to Gray Code
inline int gray_code(int x) { return x ^(x >> 1); }


/*!
  \brief Convert anything to string

  \param i (Input) The value to be converted to a string
*/
template <typename T>
std::string to_str(const T &i);

/*!
  \brief Convert double to string

  \param[in]  i          The value to be converted to a string
  \param[in]  precision  The number of digits used to represent the
  fractional part
*/
std::string to_str(const double &i, const int precision);

//@}

template <typename T>
std::string to_str(const T &i)
{
  std::ostringstream ss;
  ss.precision(8);
  ss.setf(std::ostringstream::scientific, std::ostringstream::floatfield);
  ss << i;
  return ss.str();
}

//! \cond

// ---------------------------------------------------------------------
// Instantiations
// ---------------------------------------------------------------------

#ifndef _MSC_VER

extern template bvec to_bvec(const svec &v);
extern template bvec to_bvec(const ivec &v);

extern template svec to_svec(const bvec &v);
extern template svec to_svec(const ivec &v);
extern template svec to_svec(const vec &v);

extern template ivec to_ivec(const bvec &v);
extern template ivec to_ivec(const svec &v);
extern template ivec to_ivec(const vec &v);

extern template vec to_vec(const bvec &v);
extern template vec to_vec(const svec &v);
extern template vec to_vec(const ivec &v);

extern template cvec to_cvec(const bvec &v);
extern template cvec to_cvec(const svec &v);
extern template cvec to_cvec(const ivec &v);
extern template cvec to_cvec(const vec &v);

extern template cvec to_cvec(const bvec &real, const bvec &imag);
extern template cvec to_cvec(const svec &real, const svec &imag);
extern template cvec to_cvec(const ivec &real, const ivec &imag);
extern template cvec to_cvec(const vec &real, const vec &imag);

extern template bmat to_bmat(const smat &m);
extern template bmat to_bmat(const imat &m);

extern template smat to_smat(const bmat &m);
extern template smat to_smat(const imat &m);
extern template smat to_smat(const mat &m);

extern template imat to_imat(const bmat &m);
extern template imat to_imat(const smat &m);
extern template imat to_imat(const mat &m);

extern template mat to_mat(const bmat &m);
extern template mat to_mat(const smat &m);
extern template mat to_mat(const imat &m);

extern template cmat to_cmat(const bmat &m);
extern template cmat to_cmat(const smat &m);
extern template cmat to_cmat(const imat &m);
extern template cmat to_cmat(const mat &m);

extern template cmat to_cmat(const bmat &real, const bmat &imag);
extern template cmat to_cmat(const smat &real, const smat &imag);
extern template cmat to_cmat(const imat &real, const imat &imag);
extern template cmat to_cmat(const mat &real, const mat &imag);

#endif // _MSC_VER

//! \endcond

} // namespace itpp

#endif // CONVERTERS_H
