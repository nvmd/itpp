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
  \brief Definitions of converters between different vector and matrix types
  \author Tony Ottosson, Tobias Ringström and Pål Frenger

  $Revision$

  $Date$
*/

#ifndef __converters_h
#define __converters_h

#include "itpp/base/vec.h"
#include "itpp/base/mat.h"
#include <sstream>


namespace itpp {

  /*! \addtogroup convertfunc
   */
  //!@{

  /*!
    \relates Vec
    \brief Converts a Vec<T> to bvec
  */
  template <class T>
    bvec to_bvec(const Vec<T> &v);

  /*!
    \relates Vec
    \brief Converts a Vec<T> to svec
  */
  template <class T>
    svec to_svec(const Vec<T> &v);

  /*!
    \relates Vec
    \brief Converts a Vec<T> to ivec
  */
  template <class T>
    ivec to_ivec(const Vec<T> &v);

  /*!
    \relates Vec
    \brief Converts a Vec<T> to vec
  */
  template <class T>
    vec to_vec(const Vec<T> &v);

  /*!
    \relates Vec
    \brief Converts a Vec<T> to cvec
  */
  template <class T>
    cvec to_cvec(const Vec<T> &v);

  /*!
    \relates Vec
    \brief Converts real and imaginary Vec<T> to cvec
  */
  template <class T>
    cvec to_cvec(const Vec<T> &real, const Vec<T> &imag);

  /*!
    \relates Vec
    \brief Converts an int to ivec
  */
  ivec to_ivec(int s);

  /*!
    \relates Vec
    \brief Converts an double to vec
  */
  vec to_vec(double s);

  /*!
    \relates Vec
    \brief Converts real and imaginary double to cvec
  */
  cvec to_cvec(double real, double imag);

  /*!
    \relates Mat
    \brief Converts a Mat<T> to bmat
  */
  template <class T>
    bmat to_bmat(const Mat<T> &m);

  /*!
    \relates Mat
    \brief Converts a Mat<T> to smat
  */
  template <class T>
    smat to_smat(const Mat<T> &m);

  /*!
    \relates Mat
    \brief Converts a Mat<T> to imat
  */
  template <class T>
    imat to_imat(const Mat<T> &m);

  /*!
    \relates Mat
    \brief Converts a Mat<T> to mat
  */
  template <class T>
    mat to_mat(const Mat<T> &m);

  /*!
    \relates Mat
    \brief Converts a Mat<T> to cmat
  */
  template <class T>
    cmat to_cmat(const Mat<T> &m);

  /*!
    \relates Mat
    \brief Converts real and imaginary Mat<T> to cmat
  */
  template <class T>
    cmat to_cmat(const Mat<T> &real, const Mat<T> &imag);

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

  /*! 
    \brief Convert anything to string
    
    \param i (Input) The value to be converted to a string
  */
  template <typename T>
    std::string to_str(const T &i);
  
  /*! 
    \brief Convert double to string
    
    \param i (Input) The value to be converted to a string
    \param precision (Input) The number of digits used to represent the fractional part
  */
  std::string to_str(const double &i, const int precision);

  //!@}

  template <typename T>
    std::string to_str(const T &i)
    {
      std::ostringstream ss;
	  ss.precision(8);
	  ss.setf(std::ostringstream::scientific,std::ostringstream::floatfield);
      ss << i;
      return ss.str();
    }

  // ---------------------- Instantiations -----------------------------------------
#if 0
#ifndef _MSC_VER

  //! Template instantiation of to_bvec
  extern template bvec to_bvec(const svec &v);
  //! Template instantiation of to_bvec
  //#ifndef _MSC_VER
  extern template bvec to_bvec(const Vec<int> &v);
  //#endif
  //! Template instantiation of to_svec
  extern template svec to_svec(const bvec &v);
  //! Template instantiation of to_svec
  extern template svec to_svec(const ivec &v);
  //! Template instantiation of to_svec
  extern template svec to_svec(const svec &v);

  //#ifndef _MSC_VER
  //! Template instantiation of to_ivec
  extern template ivec to_ivec(const bvec &v);
  //#endif
  //! Template instantiation of to_ivec
  extern template ivec to_ivec(const svec &v);
  //! Template instantiation of to_ivec
  extern template ivec to_ivec(const ivec &v);
  //! Template instantiation of to_ivec
  extern template ivec to_ivec(const vec &v);

  //! Template instantiation of to_vec
  extern template vec to_vec(const bvec &v);
  //! Template instantiation of to_vec
  extern template vec to_vec(const svec &v);
  //! Template instantiation of to_vec
  extern template vec to_vec(const ivec &v);
  //! Template instantiation of to_vec
  extern template vec to_vec(const vec &v);

  // Template instantiation of to_cvec
  //template cvec to_cvec(const bvec &v); //Specialization created above

  //! Template instantiation of to_cvec
  extern template cvec to_cvec(const svec &v);
  //! Template instantiation of to_cvec
  extern template cvec to_cvec(const ivec &v);
  //! Template instantiation of to_cvec
  extern template cvec to_cvec(const vec &v);
  //! Template instantiation of to_cvec
  extern template cvec to_cvec(const cvec &v);
  //! Template instantiation of to_cvec
  extern template cvec to_cvec(const bvec &real, const bvec &imag);
  //! Template instantiation of to_cvec
  extern template cvec to_cvec(const svec &real, const svec &imag);
  //! Template instantiation of to_cvec
  extern template cvec to_cvec(const ivec &real, const ivec &imag);
  //! Template instantiation of to_cvec
  extern template cvec to_cvec(const vec &real, const vec &imag);

  //! Template instantiation of to_bmat
  extern template bmat to_bmat(const smat &m);
  //! Template instantiation of to_bmat
  extern template bmat to_bmat(const imat &m);
  //! Template instantiation of to_bmat
  extern template smat to_smat(const bmat &m);
  //! Template instantiation of to_bmat
  extern template smat to_smat(const imat &m);

  //! Template instantiation of to_imat
  extern template imat to_imat(const bmat &m);
  //! Template instantiation of to_imat
  extern template imat to_imat(const smat &m);
  //! Template instantiation of to_imat
  extern template imat to_imat(const imat &m);
  // Template instantiation of to_imat
  extern template imat to_imat(const mat &m);

  //! Template instantiation of to_mat
  extern template mat to_mat(const bmat &m);
  //! Template instantiation of to_mat
  extern template mat to_mat(const smat &m);
  //! Template instantiation of to_mat
  extern template mat to_mat(const imat &m);
  //! Template instantiation of to_mat
  extern template mat to_mat(const mat &m);

  // Template instantiation of to_cmat
  //template cmat to_cmat(const bmat &m); //Specialization created above

  //! Template instantiation of to_cmat
  extern template cmat to_cmat(const smat &m);
  //! Template instantiation of to_cmat
  extern template cmat to_cmat(const imat &m);
  //! Template instantiation of to_cmat
  extern template cmat to_cmat(const mat &m);
  //! Template instantiation of to_cmat
  extern template cmat to_cmat(const cmat &m);
  //! Template instantiation of to_cmat
  extern template cmat to_cmat(const bmat &real, const bmat &imag);
  //! Template instantiation of to_cmat
  extern template cmat to_cmat(const smat &real, const smat &imag);
  //! Template instantiation of to_cmat
  extern template cmat to_cmat(const imat &real, const imat &imag);
  //! Template instantiation of to_cmat
  extern template cmat to_cmat(const mat &real, const mat &imag);
#endif
#endif

} //namespace itpp

#endif // __converters_h
