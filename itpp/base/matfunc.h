/*!
 * \file
 * \brief Definitions of functions on vectors and matrices
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

#ifndef MATFUNC_H
#define MATFUNC_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/converters.h>
#include <itpp/base/scalfunc.h>
#include <itpp/base/itassert.h>
#include <itpp/base/specmat.h>
#include <itpp/base/binary.h>
#include <itpp/base/sort.h>


namespace itpp {


  //! Define of absolute(x)
#define absolut(x) (((x)>0) ? (x) : (-(x)))

  /*!
    \addtogroup matrix_functions
    \brief Functions on vectors and matrices
  */

  //!@{

  //! Length of vector
	
	template<class T>
    int length(const Vec<T> &v) { return v.length(); }

  //! Length of vector
  template<class T>
    int size(const Vec<T> &v) { return v.length(); }



  //! Sum of all elements in the vector 
  template<class T>
    T sum(const Vec<T> &v)
    {
      T M=0;

      for (int i=0;i<v.length();i++)
	M += v[i];

      return M;
    }

  /*! 
    \brief sum of elements in the matrix \c m
    sum(m)=sum(m,1) returns a vector where the elements are sum over each column
    sum(m,2) returns a vector where the elements are sum over each row
  */
  template<class T>
    Vec<T> sum(const Mat<T> &m, int dim=1)
    {
      it_assert(dim==1 || dim==2, "sum: dimension need to be 1 or 2");
      Vec<T> out;

      if (dim == 1) {
	out.set_size(m.cols(), false);

	for (int i=0; i<m.cols(); i++)
	  out(i) = sum(m.get_col(i));
      } else {
	out.set_size(m.rows(), false);

	for (int i=0; i<m.rows(); i++)
	  out(i) = sum(m.get_row(i));
      }
      
      return out;
    }

  //! Sum of square of the elements in a vector
  template<class T>
    T sum_sqr(const Vec<T> &v)
    {
      T M=0;

      for (int i=0; i<v.length(); i++)
	M += v[i] * v[i];

      return M;
    }

  /*! 
    \brief sum of the square of elements in the matrix \c m
    sum(m)=sum(m,1) returns a vector where the elements are sum squared over each column
    sum(m,2) returns a vector where the elements are sum squared over each row
  */
  template<class T>
    Vec<T> sum_sqr(const Mat<T> &m, int dim=1)
    {
      it_assert(dim==1 || dim==2, "sum_sqr: dimension need to be 1 or 2");
      Vec<T> out;

      if (dim == 1) {
	out.set_size(m.cols(), false);

	for (int i=0; i<m.cols(); i++)
	  out(i) = sum_sqr(m.get_col(i));
      } else {
	out.set_size(m.rows(), false);

	for (int i=0; i<m.rows(); i++)
	  out(i) = sum_sqr(m.get_row(i));
      }
      
      return out;
    }

  //! Cumulative sum of all elements in the vector 
  template<class T>
    Vec<T> cumsum(const Vec<T> &v)
    {
      Vec<T> out(v.size());

      out(0)=v(0);
      for (int i=1; i<v.size(); i++)
	out(i) = out(i-1) + v(i);

      return out;
    }

  /*! 
    \brief cumulative sum of elements in the matrix \c m
    cumsum(m)=cumsum(m,1) returns a matrix where the elements are sums over each column
    cumsum(m,2) returns a matrix where the elements are sums over each row
  */
  template<class T>
    Mat<T> cumsum(const Mat<T> &m, int dim=1)
    {
      it_assert(dim==1 || dim==2, "cumsum: dimension need to be 1 or 2");
      Mat<T> out(m.rows(), m.cols());

      if (dim == 1) {
	for (int i=0; i<m.cols(); i++)
	  out.set_col(i, cumsum(m.get_col(i)));
      } else {
	for (int i=0; i<m.rows(); i++)
	  out.set_row(i, cumsum(m.get_row(i)));
      }

      return out;
    }
  
  //! The product of all elements in the vector
  template<class T>
    T prod(const Vec<T> &v)
    {
      it_assert(v.size() >= 1, "prod: size of vector should be at least 1");
      T out = v(0);

      for (int i=1; i<v.size(); i++)
	out *= v(i);

      return out;
    }

  /*! 
    \brief product of elements in the matrix \c m
    prod(m)=prod(m,1) returns a vector where the elements are products over each column
    prod(m,2) returns a vector where the elements are products over each row
  */
  template<class T>
    Vec<T> prod(const Mat<T> &m, int dim=1)
    {
      it_assert(dim==1 || dim==2, "prod: dimension need to be 1 or 2");
      Vec<T> out(m.cols());

      if (dim == 1) {
	it_assert(m.cols() >= 1 && m.rows() >= 1, "prod: number of columns should be at least 1");
	out.set_size(m.cols(), false);

	for (int i=0; i<m.cols(); i++)
	  out(i) = prod(m.get_col(i));
      } else {
	it_assert(m.cols() >= 1 && m.rows() >= 1, "prod: number of rows should be at least 1");
	out.set_size(m.rows(), false);

	for (int i=0; i<m.rows(); i++)
	  out(i) = prod(m.get_row(i));
      }
      return out;
    }

  //! Vector cross product. Vectors need to be of size 3
  template<class T>
    Vec<T> cross(const Vec<T> &v1, const Vec<T> &v2)
    {
      it_assert( v1.size() == 3 && v2.size() == 3, "cross: vectors should be of size 3");

      Vec<T> r(3);

      r(0) = v1(1) * v2(2) - v1(2) * v2(1);
      r(1) = v1(2) * v2(0) - v1(0) * v2(2);
      r(2) = v1(0) * v2(1) - v1(1) * v2(0);

      return r;
    }




  //! Apply arbitrary function to a vector
  template<class T, class fT>
    Vec<T> apply_function(fT (*f)(fT), const Vec<T> &data)
  {
    Vec<T> out(data.length());

    for (int i=0;i<data.length();i++)
      out[i]=T(f(fT(data[i])));
    return out;
  }


  //! Apply arbitrary functions to a matrix
  template<class T, class fT>
    Mat<T> apply_function(fT (*f)(fT), const Mat<T> &data)
  {
    Mat<T> out(data.rows(),data.cols());

    for (int i=0;i<out.rows();i++)
      for (int j=0;j<out.cols();j++)
	//out(i,j)=static_cast<T>(f(static_cast<fT>(data(i,j))));
	out(i,j)=T(f(fT(data(i,j))));

    return out;
  }



  //! Zero-pad a vector to size n
  template<class T>
    Vec<T> zero_pad(const Vec<T> &v, int n)
    {
      it_assert(n>=v.size(), "zero_pad() cannot shrink the vector!");
      Vec<T> v2(n);
      v2.set_subvector(0, v.size()-1, v);
      if (n > v.size())
	v2.set_subvector(v.size(), n-1, T(0));

      return v2;
    }

  //! Zero-pad a vector to the nearest greater power of two
  template<class T>
    Vec<T> zero_pad(const Vec<T> &v)
    {
      int n = pow2(needed_bits(v.size()));

      return n==v.size() ? v : zero_pad(v, n);
    }

  //! Zero-pad a matrix to size rows x cols
  template<class T>
    Mat<T> zero_pad(const Mat<T> &m, int rows, int cols)
    {
      it_assert(rows>=m.rows() && cols>=m.cols(), "zero_pad() cannot shrink the matrix!");
      Mat<T> m2(rows, cols);
      m2.set_submatrix(0,m.rows()-1,0,m.cols()-1, m);
      if (cols > m.cols()) // Zero
	m2.set_submatrix(0,m.rows()-1, m.cols(),cols-1, T(0));
      if (rows > m.rows()) // Zero
	m2.set_submatrix(m.rows(), rows-1, 0, cols-1, T(0));

      return m2;
    }
  

  //! Return zero if indexing outside the vector \c v otherwise return the element \c index
  template<class T>
  T index_zero_pad(const Vec<T> &v, const int index)
  {
    if (index >= 0 && index < v.size())
      return v(index);
    else
      return T(0);
  }

  /*! 
    \brief Transposition of the matrix \c m returning the transposed matrix in \c out
  */
  template<class T>
    void transpose(const Mat<T> &m, Mat<T> &out) { out = m.T(); }

  /*! 
    \brief Transposition of the matrix \c m
  */
  template<class T>
    Mat<T> transpose(const Mat<T> &m) { return m.T(); }

  /*! 
    \brief Hermitian transpose (complex conjugate transpose) of the matrix \c m returning the transposed matrix in \c out
  */
  template<class T>
    void hermitian_transpose(const Mat<T> &m, Mat<T> &out) { out = m.H(); }

  /*!  
    \brief Hermitian transpose (complex conjugate transpose) of the matrix \c m
  */
  template<class T>
    Mat<T> hermitian_transpose(const Mat<T> &m) { return m.H(); }


  /*! 
    \brief Computes the Kronecker product of two matrices

    <tt>k = kron(x, y)</tt> returns the Kronecker tensor product of \c
    x and \c y. The result is a large array formed by taking all
    possible products between the elements of \c x and those of \c
    Y. If \c x is <tt>(m x n)</tt> and \c y is <tt>(p x q)</tt>, then
    <tt>kron(x, y)</tt> is <tt>(m*p x n*q)</tt>.

    \author Adam Piatyszek
  */
  template<class Num_T>
    Mat<Num_T> kron(const Mat<Num_T>& x, const Mat<Num_T>& y) {
    int x_rows = x.rows();
    int x_cols = x.cols();
    int y_rows = y.rows();
    int y_cols = y.cols();
    
    Mat<Num_T> result(x_rows * y_rows, x_cols * y_cols);

    for (int i = 0; i < x_rows; i++) {
      for (int j = 0; j < x_cols; j++) {
	result.set_submatrix(i * y_rows, j * y_cols, x(i, j) * y);
      }
    }
    
    return result;
  }
  //!@}



  // -------------------- Diagonal matrix functions ---------------------------------------

  /*! 
    \addtogroup diag
    
  */

  //!@{

  /*! 
    \brief Returns a diagonal matrix whith the elements of the vector \c v on the diagonal and zeros elsewhere.

    The size of the return matrix will be \f$n \times n\f$, where \f$n\f$ is the length of the input vector \c v.
  */
  template<class T>
    Mat<T> diag(const Vec<T> &v)
    {
      Mat<T> m(v.size(), v.size());
      m = T(0);
      for (int i=v.size()-1; i>=0; i--)
	m(i,i) = v(i);
      return m;
    }

  /*! 
    \brief Returns in the output wariable \c m a diagonal matrix whith the elements of the vector \c v on the 
    diagonal and zeros elsewhere.
  
    The size of the output matrix \c m will be \f$n \times n\f$, where \f$n\f$ is the length of the input vector \c v.
  */
  template<class T>
    void diag(const Vec<T> &v, Mat<T> &m)
    {
      m.set_size(v.size(), v.size(), false);
      m = T(0);
      for (int i=v.size()-1; i>=0; i--)
	m(i,i) = v(i);
    }

  /*! 
    \brief Returns the diagonal elements of the input matrix \c m.

    The input matrix \c m must be a square \f$n \times n\f$ matrix. The size of the output vector will be \f$n\f$.
  */
  template<class T>
    Vec<T> diag(const Mat<T> &m)
    {
      Vec<T> t(std::min(m.rows(), m.cols()));

      for (int i=0; i<t.size(); i++)
	t(i) = m(i,i);

      return t;
    }

  // 
  /*!
    \brief Returns a matrix with the elements of the input vector \c main on the diagonal and the elements of 
    the input vector \c sup on the diagonal row above.

    If the number of elements in the vector \c main is \f$n\f$, then the number of elements in the input vector 
    \c sup must be \f$n-1\f$. The size of the return matrix will be \f$n \times n\f$.
  */
  template<class T>
    Mat<T> bidiag(const Vec<T> &main, const Vec<T> &sup)
    {
      it_assert(main.size() == sup.size()+1, "bidiag()");

      int n=main.size();
      Mat<T> m(n, n);
      m = T(0);
      for (int i=0; i<n-1; i++) {
	m(i,i) = main(i);
	m(i,i+1) = sup(i);
      }
      m(n-1,n-1) = main(n-1);

      return m;
    }

  /*!
    \brief Returns in the output variable \c m a matrix with the elements of the input vector \c main on the diagonal 
    and the elements of the input vector \c sup on the diagonal row above.

    If the number of elements in the vector \c main is \f$n\f$, then the number of elements in the input vector 
    \c sup must be \f$n-1\f$. The size of the output matrix \c m will be \f$n \times n\f$.
  */
  template<class T>
    void bidiag(const Vec<T> &main, const Vec<T> &sup, Mat<T> &m)
    {
      it_assert(main.size() == sup.size()+1, "bidiag()");

      int n=main.size();
      m.set_size(n, n);
      m = T(0);
      for (int i=0; i<n-1; i++) {
	m(i,i) = main(i);
	m(i,i+1) = sup(i);
      }
      m(n-1,n-1) = main(n-1);
    }

  /*!
    \brief Returns the main diagonal and the diagonal row above in the two output vectors \c main and \c sup.

    The input matrix \c in must be a square \f$n \times n\f$ matrix. The length of the output vector \c main will be \f$n\f$ 
    and the length of the output vector \c sup will be \f$n-1\f$.
  */
  template<class T>
    void bidiag(const Mat<T> &m, Vec<T> &main, Vec<T> &sup)
    {
      it_assert(m.rows() == m.cols(), "bidiag(): Matrix must be square!");

      int n=m.cols();
      main.set_size(n);
      sup.set_size(n-1);
      for (int i=0; i<n-1; i++) {
	main(i) = m(i,i);
	sup(i) = m(i,i+1);
      }
      main(n-1) = m(n-1,n-1);
    }

  /*!
    \brief Returns a matrix with the elements of \c main on the diagonal, the elements of \c sup on the diagonal row above, 
    and the elements of \c sub on the diagonal row below.

    If the length of the input vector \c main is \f$n\f$ then the lengths of the vectors \c sup and \c sub 
    must equal \f$n-1\f$. The size of the return matrix will be \f$n \times n\f$.
  */
  template<class T>
    Mat<T> tridiag(const Vec<T> &main, const Vec<T> &sup, const Vec<T> &sub)
    {
      it_assert(main.size()==sup.size()+1 && main.size()==sub.size()+1, "bidiag()");

      int n=main.size();
      Mat<T> m(n, n);
      m = T(0);
      for (int i=0; i<n-1; i++) {
	m(i,i) = main(i);
	m(i,i+1) = sup(i);
	m(i+1,i) = sub(i);
      }
      m(n-1,n-1) = main(n-1);

      return m;
    }

  /*!
    \brief Returns in the output matrix \c m a matrix with the elements of \c main on the diagonal, the elements of \c sup on the
    diagonal row above, and the elements of \c sub on the diagonal row below.
  
    If the length of the input vector \c main is \f$n\f$ then the lengths of the vectors \c sup and \c sub 
    must equal \f$n-1\f$. The size of the output matrix \c m will be \f$n \times n\f$.
  */
  template<class T>
    void tridiag(const Vec<T> &main, const Vec<T> &sup, const Vec<T> &sub, Mat<T> &m)
    {
      it_assert(main.size()==sup.size()+1 && main.size()==sub.size()+1, "bidiag()");

      int n=main.size();
      m.set_size(n, n);
      m = T(0);
      for (int i=0; i<n-1; i++) {
	m(i,i) = main(i);
	m(i,i+1) = sup(i);
	m(i+1,i) = sub(i);
      }
      m(n-1,n-1) = main(n-1);
    }

  /*!
    \brief Returns the main diagonal, the diagonal row above, and the diagonal row below int the output vectors \c main, \c sup, and \c sub.

    The input matrix \c m must be a square \f$n \times n\f$ matrix. The length of the output vector \c main will be \f$n\f$ 
    and the length of the output vectors \c sup and \c sup will be \f$n-1\f$.
  */
  template<class T>
    void tridiag(const Mat<T> &m, Vec<T> &main, Vec<T> &sup, Vec<T> &sub)
    {
      it_assert(m.rows() == m.cols(), "tridiag(): Matrix must be square!");

      int n=m.cols();
      main.set_size(n);
      sup.set_size(n-1);
      sub.set_size(n-1);
      for (int i=0; i<n-1; i++) {
	main(i) = m(i,i);
	sup(i) = m(i,i+1);
	sub(i) = m(i+1,i);
      }
      main(n-1) = m(n-1,n-1);
    }



  /*! 
    \brief The trace of the matrix \c m, i.e. the sum of the diagonal elements.
  */
  template<class T>
    T trace(const Mat<T> &m)
    {
      return sum(diag(m));
    }

  //!@}


  // ----------------- reshaping vectors and matrices ---------------------------
  /*! 
    \addtogroup reshaping
    
  */

  //!@{

  //! Reverse the input vector
  template<class T>
    Vec<T> reverse(const Vec<T> &in)
    {
      int i, s=in.length();

      Vec<T> out(s);
      for (i=0;i<s;i++)
	out[i]=in[s-1-i];
      return out;
    }

  //! Row vectorize the matrix [(0,0) (0,1) ... (N-1,N-2) (N-1,N-1)]
  template<class T>
    Vec<T> rvectorize(const Mat<T> &m)
    {
      int i, j, n=0, r=m.rows(), c=m.cols();
      Vec<T> v(r * c);

      for (i=0; i<r; i++)
	for (j=0; j<c; j++)
	  v(n++) = m(i,j);

      return v;
    }

  //! Column vectorize the matrix [(0,0) (1,0) ... (N-2,N-1) (N-1,N-1)]
  template<class T>
    Vec<T> cvectorize(const Mat<T> &m)
    {
      int i, j, n=0, r=m.rows(), c=m.cols();
      Vec<T> v(r * c);

      for (j=0; j<c; j++)
	for (i=0; i<r; i++)
	  v(n++) = m(i,j);

      return v;
    }

  /*!
    \brief Reshape the matrix into an rows*cols matrix

    The data is taken columnwise from the original matrix and written columnwise into the new matrix.
  */
  template<class T>
    Mat<T> reshape(const Mat<T> &m, int rows, int cols)
    {
      it_assert1(m.rows()*m.cols() == rows*cols, "Mat<T>::reshape: Sizes must match");
      Mat<T> temp(rows, cols);
      int i, j, ii=0, jj=0;
      for (j=0; j<m.cols(); j++) {
	for (i=0; i<m.rows(); i++) {
	  temp(ii++,jj) = m(i,j);
	  if (ii == rows) {
	    jj++; ii=0;
	  }
	}
      }
      return temp;
    }

  /*!
    \brief Reshape the vector into an rows*cols matrix

    The data is element by element from the vector and written columnwise into the new matrix.
  */
  template<class T>
    Mat<T> reshape(const Vec<T> &v, int rows, int cols)
    {
      it_assert1(v.size() == rows*cols, "Mat<T>::reshape: Sizes must match");
      Mat<T> temp(rows, cols);
      int i, j, ii=0;
      for (j=0; j<cols; j++) {
	for (i=0; i<rows; i++) {
	  temp(i,j) = v(ii++);
	}
      }
      return temp;
    }

  //!@}


  /*! 
    \addtogroup upsample
  */

  //!@{

  //! Repeat each element in the vector norepeats times in sequence
  template<class T>
    Vec<T> repeat(const Vec<T> &v, int norepeats)
    {
      Vec<T> temp(v.length()*norepeats);

      for(int i=0; i<v.length(); i++) {
	for(int j=0;j<norepeats;j++)
	  temp(i*norepeats+j)=v(i);
      }
      return temp;
    }

  //! Repeats each column norepeats times in sequence
  template<class T>
    Mat<T> repeat(const Mat<T> &m, int norepeats)
    {
      Mat<T> temp(m.rows(), m.cols()*norepeats);

      for (int j=0; j<m.cols(); j++) {
	for (int i=0;i<norepeats;i++) {
	  temp.set_col(j*norepeats+i, m.get_col(j));
	}
      }
      return temp;
    }

  //! Upsample a vector by inserting \a (usf-1) zeros after each sample
  template<class T> 
    void upsample(const Vec<T> &v, int usf, Vec<T> &u)
    {
      it_assert1(usf >= 1, "upsample: upsampling factor must be equal or greater than one" );
      u.set_size(v.length()*usf);
      u.clear();
      for(long i=0;i<v.length();i++)
	u(i*usf)=v(i);
    }


  //! Upsample a vector by incerting \a (usf-1) zeros after each sample
  template<class T>
    Vec<T> upsample(const Vec<T> &v, int usf)
    {
      Vec<T> u;
      upsample(v,usf,u);
      return u;
    }

  //! Upsample each column by incerting \a (usf-1) zeros after each column
  template<class T>
    void upsample(const Mat<T> &v, int usf, Mat<T> &u)
    {
      it_assert1(usf >= 1, "upsample: upsampling factor must be equal or greater than one" );
      u.set_size(v.rows(),v.cols()*usf);
      u.clear();
      for (long j=0;j<v.cols();j++)
	u.set_col(j*usf,v.get_col(j));
    }

  //! Upsample each column by incerting \a (usf-1) zeros after each column
  template<class T>
    Mat<T> upsample(const Mat<T> &v, int usf)
    {
      Mat<T> u;
      upsample(v,usf,u);
      return u;
    }

  //! Upsample each column by a factor of  \a (usf-1) by linear interpolation
  template<class T>
    void lininterp(const Mat<T> &m, int usf, Mat<T> &u)
    {
      it_assert1(usf >= 1, "lininterp: upsampling factor must be equal or greater than one" );
      long L = (m.cols()-1)*usf+1;
      u.set_size(m.rows(),L);
      for (long i = 0; i < m.rows(); i++){
	for (long j = 0; j < L-1; j++)
	  u(i,j) = (m(i,j/usf) + (j % usf)/((double)usf)*(m(i,(j+usf)/usf)-m(i,j/usf)));
	u(i,L-1) = m(i,m.cols()-1);
      }
    }

  /*! 
    \brief Upsample each column of matrix \a m to achieve \a f_ups
    frequency using linear interpolation

    This function performs upsampling of matrix \a m to achieve \a
    nrof_samples samples at \a f_ups frequency starting from the 
    sample at \a t_start time. The frequency of input samples stored
    in the matrix \a m is defined by the \a f_base parameter.

    \author Adam Piatyszek
  */
  template<class T>
    Mat<T> lininterp(const Mat<T> &m, const double f_base, const double f_ups, const int nrof_samples, const double t_start = 0)
    {
      double t_base = 1 / f_base;
      double t_ups = 1 / f_ups;
      int rows = m.rows();
      int cols = m.cols();
      it_assert1(f_ups > f_base, "lininterp: upsampled frequency must be greater than base frequency" );
      it_assert1((t_start >= 0) && (t_start < cols * t_base), "lininterp: incorrect start time offset");
      it_assert1((nrof_samples * t_ups + t_start) <= (cols * t_base), "lininterp: too many samples required or input data to short");
      Mat<T> u(rows, nrof_samples);
      double curr_time = t_start;
    
      int i = 0;
      int k = 0;
      while (i < cols - 1) {
	while ((curr_time < (i + 1) * t_base) && (k < nrof_samples)) {
	  for (int j = 0; j < rows; j++) {
	    u(j, k) = (m(j, i) * ((i + 1) * t_base - curr_time) 
		       - m(j, i + 1) * (i * t_base - curr_time)) / t_base;
	  }
	  k++;
	  curr_time += t_ups;
	}
	i++;
      }
      return u;
    }


  //! Upsample each column by a factor of  \a (usf-1) by linear interpolation
  template<class T>
    Mat<T> lininterp(const Mat<T> &m, int usf)
    {
      Mat<T> u;
      lininterp(m,usf,u);
      return u;
    }

  //! Upsample by a factor of  \a (usf-1) by linear interpolation
  template<class T>
    void lininterp(const Vec<T> &v, int usf, Vec<T> &u)
    {
      it_assert1(usf >= 1, "lininterp: upsampling factor must be equal or greater than one" );
      long L = (v.length()-1)*usf+1;
      u.set_size(L);
      for (long j = 0; j < L-1; j++) {
	u(j) = (v(j/usf) + (j % usf)/((double)usf)*(v((j+usf)/usf)-v(j/usf)));
      }
      u(L-1) = v(v.length()-1);
    }

  //! Upsample by a factor of  \a (usf-1) by linear interpolation
  template<class T>
    Vec<T> lininterp(const Vec<T> &v, int usf)
    {
      Vec<T> u;
      lininterp(v,usf,u);
      return u;
    }

  /*! 
    \brief Upsample each sample of vector \a v to achieve \a f_ups
    frequency using linear interpolation

    This function performs upsampling of vector \a v to achieve \a
    nrof_samples samples at \a f_ups frequency starting from the 
    sample at \a t_start time. The frequency of input samples stored
    in the vector \a v is defined by the \a f_base parameter.

    \author Adam Piatyszek
  */
  template<class T>
    Vec<T> lininterp(const Vec<T> &v, const double f_base, const double f_ups, const int nrof_samples, const double t_start = 0)
    {
      double t_base = 1 / f_base;
      double t_ups = 1 / f_ups;
      int len = v.length();
      it_assert1(f_ups > f_base, "lininterp: upsampled frequency must be greater than base frequency" );
      it_assert1((t_start >= 0) && (t_start < len * t_base), "lininterp: incorrect start time offset");
      it_assert1((nrof_samples * t_ups + t_start) <= (len * t_base), "lininterp: too many samples required or input data to short");
      Vec<T> u(nrof_samples);
      double curr_time = t_start;

      int i = 0;
      int k = 0;
      while (i < len - 1) {
	while ((curr_time < (i + 1) * t_base) && (k < nrof_samples)) {
	  u(k) = (v(i) * ((i + 1) * t_base - curr_time) 
		  - v(i + 1) * (i * t_base - curr_time)) / t_base;
	  k++;
	  curr_time += t_ups;
	}
	i++;
      }
      return u;
    }

  //!@}

  // ---------------------- Instantiations -----------------------------------------
#ifndef _MSC_VER

  //! Extern Template instantiation of length
  extern template int length(const vec &v);
  //! Extern Template instantiation of length
  extern template int length(const cvec &v);
  //! Extern Template instantiation of length
  extern template int length(const svec &v);
  //! Extern Template instantiation of length
  extern template int length(const ivec &v);
  //! Extern Template instantiation of length
  extern template int length(const bvec &v);

  //! Extern Template instantiation of sum
  extern template double sum(const vec &v);
  //! Extern Template instantiation of sum
  extern template std::complex<double> sum(const cvec &v);
  //! Extern Template instantiation of sum
  extern template short sum(const svec &v);
  //! Extern Template instantiation of sum
  extern template int sum(const ivec &v);
  //! Extern Template instantiation of sum
  extern template bin sum(const bvec &v);

  //! Extern Template instantiation of sum_sqr
  extern template double sum_sqr(const vec &v);
  //! Extern Template instantiation of sum_sqr
  extern template std::complex<double> sum_sqr(const cvec &v);
  //! Extern Template instantiation of sum_sqr
  extern template short sum_sqr(const svec &v);
  //! Extern Template instantiation of sum_sqr
  extern template int sum_sqr(const ivec &v);
  //! Extern Template instantiation of sum_sqr
  extern template bin sum_sqr(const bvec &v);

  //! Extern Template instantiation of cumsum
  extern template vec cumsum(const vec &v);
  //! Extern Template instantiation of cumsum
  extern template cvec cumsum(const cvec &v);
  //! Extern Template instantiation of cumsum
  extern template svec cumsum(const svec &v);
  //! Extern Template instantiation of cumsum
  extern template ivec cumsum(const ivec &v);
  //! Extern Template instantiation of cumsum
  extern template bvec cumsum(const bvec &v);

  //! Extern Template instantiation of product
  extern template double prod(const vec &v);
  //! Extern Template instantiation of product
  extern template std::complex<double> prod(const cvec &v);
  //! Extern Template instantiation of product
  extern template short prod(const svec &v);
  //! Extern Template instantiation of product
  extern template int prod(const ivec &v);
  //! Extern Template instantiation of product
  extern template bin prod(const bvec &v);

  //! Extern Template instantiation of cross
  extern template vec cross(const vec &v1, const vec &v2);
  //! Extern Template instantiation of cross
  extern template ivec cross(const ivec &v1, const ivec &v2);
  //! Extern Template instantiation of cross
  extern template svec cross(const svec &v1, const svec &v2);

  //! Extern Template instantiation of reverse
  extern template vec reverse(const vec &in);
  //! Extern Template instantiation of reverse
  extern template cvec reverse(const cvec &in);
  //! Extern Template instantiation of reverse
  extern template svec reverse(const svec &in);
  //! Extern Template instantiation of reverse
  extern template ivec reverse(const ivec &in);
  //! Extern Template instantiation of reverse
  extern template bvec reverse(const bvec &in);

  //! Extern Template instantiation of repeat
  extern template vec repeat(const vec &v, int norepeats);
  //! Extern Template instantiation of repeat
  extern template cvec repeat(const cvec &v, int norepeats);
  //! Extern Template instantiation of repeat
  extern template svec repeat(const svec &v, int norepeats);
  //! Extern Template instantiation of repeat
  extern template ivec repeat(const ivec &v, int norepeats);
  //! Extern Template instantiation of repeat
  extern template bvec repeat(const bvec &v, int norepeats);

  //! Extern Template instantiation of apply_function
  extern template vec apply_function(float (*f)(float), const vec &data);
  //! Extern Template instantiation of apply_function
  extern template vec apply_function(double (*f)(double), const vec &data);
  //! Extern Template instantiation of apply_function
  extern template cvec apply_function(std::complex<double> (*f)(std::complex<double>), const cvec &data);
  //! Extern Template instantiation of apply_function
  extern template svec apply_function(short (*f)(short), const svec &data);
  //! Extern Template instantiation of apply_function
  extern template ivec apply_function(int (*f)(int), const ivec &data);
  //! Extern Template instantiation of apply_function
  extern template bvec apply_function(bin (*f)(bin), const bvec &data);

  //! Extern Template instantiation of zero_pad
  extern template ivec zero_pad(const ivec &v, int n);
  //! Extern Template instantiation of zero_pad
  extern template vec zero_pad(const vec &v, int n);
  //! Extern Template instantiation of zero_pad
  extern template cvec zero_pad(const cvec &v, int n);
  //! Extern Template instantiation of zero_pad
  extern template bvec zero_pad(const bvec &v, int n);

  //! Extern Template instantiation of zero_pad
  extern template ivec zero_pad(const ivec &v);
  //! Extern Template instantiation of zero_pad
  extern template vec zero_pad(const vec &v);
  //! Extern Template instantiation of zero_pad
  extern template cvec zero_pad(const cvec &v);
  //! Extern Template instantiation of zero_pad
  extern template bvec zero_pad(const bvec &v);

  //! Extern Template instantiation of zero_pad
  extern template mat  zero_pad(const mat &, int, int);
  //! Extern Template instantiation of zero_pad
  extern template cmat zero_pad(const cmat &, int, int);
  //! Extern Template instantiation of zero_pad
  extern template imat zero_pad(const imat &, int, int);
  //! Extern Template instantiation of zero_pad
  extern template bmat zero_pad(const bmat &, int, int);

  //! Extern Template instantiation of sum
  extern template vec sum(const mat &m, int dim);
  //! Extern Template instantiation of sum
  extern template cvec sum(const cmat &m, int dim);
  //! Extern Template instantiation of sum
  extern template svec sum(const smat &m, int dim);
  //! Extern Template instantiation of sum
  extern template ivec sum(const imat &m, int dim);
  //! Extern Template instantiation of sum
  extern template bvec sum(const bmat &m, int dim);

  //! Extern Template instantiation of sum_sqr
  extern template vec sum_sqr(const mat & m, int dim);
  //! Extern Template instantiation of sum_sqr
  extern template cvec sum_sqr(const cmat &m, int dim);
  //! Extern Template instantiation of sum_sqr
  extern template svec sum_sqr(const smat &m, int dim);
  //! Extern Template instantiation of sum_sqr
  extern template ivec sum_sqr(const imat &m, int dim);
  //! Extern Template instantiation of sum_sqr
  extern template bvec sum_sqr(const bmat &m, int dim);

  //! Extern Template instantiation of cumsum
  extern template mat cumsum(const mat &m, int dim);
  //! Extern Template instantiation of cumsum
  extern template cmat cumsum(const cmat &m, int dim);
  //! Extern Template instantiation of cumsum
  extern template smat cumsum(const smat &m, int dim);
  //! Extern Template instantiation of cumsum
  extern template imat cumsum(const imat &m, int dim);
  //! Extern Template instantiation of cumsum
  extern template bmat cumsum(const bmat &m, int dim);

  //! Extern Template instantiation of product
  extern template vec prod(const mat &m, int dim);
  // Extern Template instantiation of product
  extern template cvec prod(const cmat &v, int dim);
  //! Extern Template instantiation of product
  extern template svec prod(const smat &m, int dim);
  //! Extern Template instantiation of product
  extern template ivec prod(const imat &m, int dim);

  //! Extern Template instantiation of diag
  extern template vec diag(const mat &in);
  //! Extern Template instantiation of diag
  extern template cvec diag(const cmat &in);

  //! Extern Template instantiation of diag
  extern template void diag(const vec &in, mat &m);
  //! Extern Template instantiation of diag
  extern template void diag(const cvec &in, cmat &m);

  //! Extern Template instantiation of diag
  extern template mat diag(const vec &v);
  //! Extern Template instantiation of diag
  extern template cmat diag(const cvec &v);

  //! Extern Template instantiation of bidiag
  extern template mat bidiag(const vec &, const vec &);
  //! Extern Template instantiation of bidiag
  extern template cmat bidiag(const cvec &, const cvec &);

  //! Extern Template instantiation of bidiag
  extern template void bidiag(const vec &, const vec &, mat &);
  //! Extern Template instantiation of bidiag
  extern template void bidiag(const cvec &, const cvec &, cmat &);

  //! Extern Template instantiation of bidiag
  extern template void bidiag(const mat &, vec &, vec &);
  //! Extern Template instantiation of bidiag
  extern template void bidiag(const cmat &, cvec &, cvec &);

  //! Extern Template instantiation of tridiag
  extern template mat tridiag(const vec &main, const vec &, const vec &);
  //! Extern Template instantiation of tridiag
  extern template cmat tridiag(const cvec &main, const cvec &, const cvec &);

  //! Extern Template instantiation of tridiag
  extern template void tridiag(const vec &main, const vec &, const vec &, mat &);
  //! Extern Template instantiation of tridiag
  extern template void tridiag(const cvec &main, const cvec &, const cvec &, cmat &);

  //! Extern Template instantiation of tridiag
  extern template void tridiag(const mat &m, vec &, vec &, vec &);
  //! Extern Template instantiation of tridiag
  extern template void tridiag(const cmat &m, cvec &, cvec &, cvec &);

  //! Extern Template instantiation of trace
  extern template double trace(const mat &in);
  //! Extern Template instantiation of trace
  extern template std::complex<double> trace(const cmat &in);
  //! Extern Template instantiation of trace
  extern template short trace(const smat &in);
  //! Extern Template instantiation of trace
  extern template int trace(const imat &in);
  //! Extern Template instantiation of trace
  extern template bin trace(const bmat &in);

  //! Extern Template instantiation of transpose
  extern template void transpose(const mat &m, mat &out);
  //! Extern Template instantiation of transpose
  extern template void transpose(const cmat &m, cmat &out);
  //! Extern Template instantiation of transpose
  extern template void transpose(const smat &m, smat &out);
  //! Extern Template instantiation of transpose
  extern template void transpose(const imat &m, imat &out);
  //! Extern Template instantiation of transpose
  extern template void transpose(const bmat &m, bmat &out);

  //! Extern Template instantiation of transpose
  extern template mat transpose(const mat &m);
  //! Extern Template instantiation of transpose
  extern template cmat transpose(const cmat &m);
  //! Extern Template instantiation of transpose
  extern template smat transpose(const smat &m);
  //! Extern Template instantiation of transpose
  extern template imat transpose(const imat &m);
  //! Extern Template instantiation of transpose
  extern template bmat transpose(const bmat &m);


  //! Extern Template instantiation of hermitian transpose
  extern template void hermitian_transpose(const mat &m, mat &out);
  //! Extern Template instantiation of hermitian transpose
  extern template void hermitian_transpose(const cmat &m, cmat &out);
  //! Extern Template instantiation of hermitian transpose
  extern template void hermitian_transpose(const smat &m, smat &out);
  //! Extern Template instantiation of hermitian transpose
  extern template void hermitian_transpose(const imat &m, imat &out);
  //! Extern Template instantiation of hermitian transpose
  extern template void hermitian_transpose(const bmat &m, bmat &out);

  //! Extern Template instantiation of hermitian transpose
  extern template mat hermitian_transpose(const mat &m);
  //! Extern Template instantiation of hermitian transpose
  extern template cmat hermitian_transpose(const cmat &m);
  //! Extern Template instantiation of hermitian transpose
  extern template smat hermitian_transpose(const smat &m);
  //! Extern Template instantiation of hermitian transpose
  extern template imat hermitian_transpose(const imat &m);
  //! Extern Template instantiation of hermitian transpose
  extern template bmat hermitian_transpose(const bmat &m);

  //! Extern Template instantiation of repeat
  extern template mat repeat(const mat &m, int norepeats);
  //! Extern Template instantiation of repeat
  extern template cmat repeat(const cmat &m, int norepeats);
  //! Extern Template instantiation of repeat
  extern template smat repeat(const smat &m, int norepeats);
  //! Extern Template instantiation of repeat
  extern template imat repeat(const imat &m, int norepeats);
  //! Extern Template instantiation of repeat
  extern template bmat repeat(const bmat &m, int norepeats);

  //! Extern Template instantiation of apply_function
  extern template mat apply_function(float (*f)(float), const mat &data);
  //! Extern Template instantiation of apply_function
  extern template mat apply_function(double (*f)(double), const mat &data);
  //! Extern Template instantiation of apply_function
  extern template cmat apply_function(std::complex<double> (*f)(std::complex<double>), const cmat &data);
  //! Extern Template instantiation of apply_function
  extern template smat apply_function(short (*f)(short), const smat &data);
  //! Extern Template instantiation of apply_function
  extern template imat apply_function(int (*f)(int), const imat &data);
  //! Extern Template instantiation of apply_function
  extern template bmat apply_function(bin (*f)(bin), const bmat &data);

  //! Extern Template instantiation of rvectorize
  extern template  vec rvectorize(const  mat &m);
  //! Extern Template instantiation of rvectorize
  extern template cvec rvectorize(const cmat &m);
  //! Extern Template instantiation of rvectorize
  extern template  ivec rvectorize(const  imat &m);
  //! Extern Template instantiation of rvectorize
  extern template  bvec rvectorize(const  bmat &m);

  //! Extern Template instantiation of cvectorize
  extern template  vec cvectorize(const  mat &m);
  //! Extern Template instantiation of cvectorize
  extern template cvec cvectorize(const cmat &m);
  //! Extern Template instantiation of cvectorize
  extern template  ivec cvectorize(const  imat &m);
  //! Extern Template instantiation of cvectorize
  extern template  bvec cvectorize(const  bmat &m);

  //! Extern Template instantiation of reshape
  extern template  mat reshape(const  mat &m, int rows, int cols);
  //! Extern Template instantiation of reshape
  extern template cmat reshape(const cmat &m, int rows, int cols);
  //! Extern Template instantiation of reshape
  extern template  imat reshape(const  imat &m, int rows, int cols);
  //! Extern Template instantiation of reshape
  extern template  bmat reshape(const  bmat &m, int rows, int cols);

  //! Extern Template instantiation of reshape
  extern template  mat reshape(const  vec &m, int rows, int cols);
  //! Extern Template instantiation of reshape
  extern template cmat reshape(const cvec &m, int rows, int cols);
  //! Extern Template instantiation of reshape
  extern template  imat reshape(const  ivec &m, int rows, int cols);
  //! Extern Template instantiation of reshape
  extern template  bmat reshape(const  bvec &m, int rows, int cols);

  //! Extern Template instantiation of upsample
  extern template vec upsample(const vec &v, int usf);
  //! Extern Template instantiation of upsample
  extern template cvec upsample(const cvec &v, int usf);
  //! Extern Template instantiation of upsample
  extern template svec upsample(const svec &v, int usf);
  //! Extern Template instantiation of upsample
  extern template ivec upsample(const ivec &v, int usf);
  //! Extern Template instantiation of upsample
  extern template bvec upsample(const bvec &v, int usf);

  //! Extern Template instantiation of upsample
  extern template mat upsample(const mat &v, int usf);
  //! Extern Template instantiation of upsample
  extern template cmat upsample(const cmat &v, int usf);
  //! Extern Template instantiation of upsample
  extern template smat upsample(const smat &v, int usf);
  //! Extern Template instantiation of upsample
  extern template imat upsample(const imat &v, int usf);
  //! Extern Template instantiation of upsample
  extern template bmat upsample(const bmat &v, int usf);

  //! Extern Template instantiation of upsample
  extern template void upsample(const vec &v, int usf,  vec &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const cvec &v, int usf,  cvec &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const svec &v, int usf,  svec &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const ivec &v, int usf,  ivec &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const bvec &v, int usf,  bvec &u);

  //! Extern Template instantiation of upsample
  extern template void upsample(const mat &v, int usf,  mat &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const cmat &v, int usf,  cmat &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const smat &v, int usf,  smat &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const imat &v, int usf,  imat &u);
  //! Extern Template instantiation of upsample
  extern template void upsample(const bmat &v, int usf,  bmat &u);

  //! Extern Template instantiation of liniterp
  extern template vec lininterp(const vec &v, int usf);
  //! Extern Template instantiation of liniterp
  extern template cvec lininterp(const cvec &v, int usf);

  //! Extern Template instantiation of liniterp
  extern template mat lininterp(const mat &v, int usf);
  //! Extern Template instantiation of liniterp
  extern template cmat lininterp(const cmat &v, int usf);

  //! Extern Template instantiation of liniterp
  extern template void lininterp(const vec &v, int usf,  vec &u);
  //! Extern Template instantiation of liniterp
  extern template void lininterp(const cvec &v, int usf,  cvec &u);

  //! Extern Template instantiation of liniterp
  extern template void lininterp(const mat &v, int usf,  mat &u);
  //! Extern Template instantiation of liniterp
  extern template void lininterp(const cmat &v, int usf,  cmat &u);

  //! Extern Template instantiation of liniterp
  extern template mat lininterp(const mat &m, const double f_base, const double f_ups, const int nrof_samples, const double t_start);
  //! Extern Template instantiation of liniterp
  extern template cmat lininterp(const cmat &m, const double f_base, const double f_ups, const int nrof_samples, const double t_start);

  //! Extern Template instantiation of liniterp
  extern template vec lininterp(const vec &v, const double f_base, const double f_ups, const int nrof_samples, const double t_start);
  //! Extern Template instantiation of liniterp
  extern template cvec lininterp(const cvec &v, const double f_base, const double f_ups, const int nrof_samples, const double t_start);
#endif

} // namespace itpp

#endif // #ifndef MATFUNC_H

