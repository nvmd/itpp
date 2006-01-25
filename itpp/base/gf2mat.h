/*!
 * \file
 * \brief Definition of a class for algebra on GF(2) (binary) matrices
 * \author Erik G. Larsson
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
 *
 * Two representations are offered: GF2mat_sparse for sparse GF(2)
 * matrices and GF2mat for dense GF(2) matrices. Conversions between
 * dense and sparse GF(2) are also possible.
 *
 * Binary vectors are represented either via the bvec class (memory
 * typically is not an issue here) or as n*1 (or 1*n) GF(2) matrix.
 *
 * Note that the \c bmat class also provides some functionality for
 * matrix algebra over GF(2) but this class is based on \c Mat<> which
 * has a fundamentally different addressing mechanism and which is
 * much less memory efficient (\c Mat<> uses one byte memory minimum
 * per element).
 */

/*
  TODO:

  - There is probably potential for improving the efficiency of some
  of the matrix/vector multiplication operators.  (See specific
  comments in the code.)

  - Make the data representation dynamic to account for the
  possibility (?) that a short integer can store more than 16 bits.
  Or use "int" or "long int"?

  - Implement a general GF(2)-matrix class which offers the user a
  transparent interface (user does not need to know whether the
  sparse or the dense representation is used)?
*/

#ifndef GF2MAT_H
#define GF2MAT_H

#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/svec.h>
#include <itpp/base/smat.h>
#include <itpp/base/itfile.h>


/* 
   An "unsigned short int" is used to represent the data
   structure. Note that, an unsigned short int has at least 16 bits,
   possibly more depending on the implementation.  Here 16 bits is
   assumed.

   Note: the data structure is implemented as unsigned int because
   some operations rely on right shifting (for signed integers the
   right shift operator is implementation defined).
*/

// log2(number of bits in one "short int"). This is used to perform
// division. 
#define lImax 4
// (1 << (lImax - 1)). This is used as mask for (1 << lImax), when
// computing remainder  
#define mImax 15

namespace itpp {

  //! Sparse GF(2) vector
  typedef Sparse_Vec<bin> GF2vec_sparse;

  //! Sparse GF(2) matrix
  typedef Sparse_Mat<bin> GF2mat_sparse;

  /*! 
    \brief Class for dense GF(2) matrices

    This class is used to represent relatively large GF(2) matrices
    making efficient use of computer memory. (One bit in the matrix
    requires one bit of memory.)

    The \c GF2mat class extends the \c bmat class by offering a more
    memory efficient representation of data and providing several
    functios for linear algebra.  \c GF2mat may be used in lieu of \c
    bmat whenever the linear algebra functionality or improved memory
    efficiency is desired.

    See also \c GF2mat_sparse which offers an efficient representation
    of sparse GF(2) matrices.

    \relates bmat
  */
  class GF2mat {
  public:
  
    // ----------- Constructors -----------

    //! Default constructor (this gives a 1*1 matrix consisting of a single zero)
    GF2mat();
  
    //! Construct an empty (all-zero) m*n matrix 
    GF2mat(int m, int n);

    //! Construct a dense GF(2) matrix from a sparse GF(2) matrix
    GF2mat(const GF2mat_sparse &X);
  
    //! Construct a dense GF(2) matrix from a subset (m1,n1) to (m2,n2) of sparse GF(2) matrix 
    GF2mat(const GF2mat_sparse &X, int m1, int n1, int m2, int n2);  

    //! Construct a dense GF(2) matrix from a subset of columns in sparse GF(2) matrix 
    GF2mat(const GF2mat_sparse &X, ivec &columns);

    /*!  
      \brief Create a dense GF(2) matrix from a single vector. 

      \param x The input vector

      \param rc A parameter that  indicates whether the result should be a row vector (rc=0),
      or a column vector (rc=1, default)
    */
    GF2mat(const bvec &x, bool rc=1);

    //! Create a dense GF(2) matrix from a bmat
    GF2mat(const bmat &X);

    //! Create a sparse GF(2) matrix from a dense GF(2) matrix
    GF2mat_sparse sparsify();

    //! Create a bvec from a GF(2) matrix (must have one column or one row)
    bvec bvecify();

    // ----------- Elementwise manipulation and simple functions -------------

    //! Getting element
    inline bin get(int i, int j) const;

    //! Getting element
    inline bin operator()(int i, int j) const { return get(i,j); };

    //! Set element i,j to s (0 or 1)
    inline void set(int i, int j, bin s);
  
    //! Add s (0 or 1) to element (i,j)
    inline void addto_element(int i, int j, bin s);

    //! Set column j to a binary vector x
    void set_col(int j, bvec x);

    //! Set row i to a binary vector x
    void set_row(int i, bvec x);

    //! Check whether the matrix is identical to zero
    bool is_zero();

    //! Swap rows i and j
    void swap_rows(int i, int j);

    //! Swap columns i and j
    void swap_cols(int i, int j);

    /*! \brief Multiply from left with permutation matrix (permute rows).

    \param perm Permutation vector

    \param I Parameter that determines permutation.   I=0: apply permutation, I=1: apply inverse permutation
    */
    void permute_rows(ivec &perm, bool I);

    /*! \brief Multiply a matrix from right with a permutation matrix (i.e.,
      permute the columns).  

    \param I Parameter that determines permutation.   I=0: apply permutation, I=1: apply inverse permutation
    */
    void permute_cols(ivec &perm, bool I); 

    //! Transpose 
    GF2mat transpose() const;

    //! Submatrix from (m1,n1) to (m2,n2)
    GF2mat get_submatrix(int m1, int n1, int m2, int n2);

    //! Concatenate horizontally (append X on the right side of matrix)
    GF2mat concatenate_horizontal(const GF2mat &X);

    //! Concatenate vertically (append X underneath)
    GF2mat concatenate_vertical(const GF2mat &X);
  
    //! Get row
    bvec get_row(int i);
  
    //! Get column
    bvec get_col(int j);
  
    //! Compute the matrix density (fraction of elements equal to "1")
    double density() const;

    //! Get number of rows
    inline int rows() { return nrows; }

    //! Get number of columns
    inline int cols() { return ncols; }

    //! Add (or equivalently, subtract) row j to row i
    void add_rows(int i, int j);

    // ---------- Linear algebra --------------

    /*!  \brief Inversion. 

    The matrix must be invertible, otherwise the
      function will terminate with an error.  
    */ 
    GF2mat inverse();

    //! Returns the number of linearly independent rows 
    int row_rank();

    /*! \brief TXP factorization. 

    Given X, compute a factorization of the form U=TXP, where U is
      upper triangular, T is square and invertible, and P is a
      permutation matrix.  This is basically an "LU"-factorization,
      but not called so here because T is not necessarily lower
      trianglular.  The matrix X may have any dimension. The
      permutation matrix P is represented as a permutation vector.

      The function returns the row rank of X.  (If X is full row rank,
      then the number of rows is returned.)
    */
    int T_fact(GF2mat &T, GF2mat &U, ivec &P);
  
    /*! \brief TXP factorization update, when bit is flipped.

      Update upper triangular factor U in the T-factorization (U=TXP)
      when the bit at position (r,c) is changed (0->1 or 1->0). The
      purpose of this function is to avoid re-running a complete
      T-factorization when a single bit is changed. The function
      assumes that T, U, and P already exist and that U=TXP before the
      function is called. The function also assumes that the rank
      provided is correct. The function updates T, U and P these
      matrices.  The function returns the new rank of the matrix after
      the bitflip.

      Note: T, U, P and the rank value supplied to the function must
      be correct one. This is not checked by the function (for reasons
      of efficiency).

      This function was primarily written specifically for use by
      certain optimization algorithms in the LDPC codec class.
    */
    int T_fact_update_bitflip(GF2mat &T, GF2mat &U, ivec &P, int rank, int r, int c);

    /*!  \brief TXP factorization update, when column is added

    Update upper triangular factor U in the T-factorization (U=TXP)
      when a column (newcol) is appended at the right side of the
      matrix.  The purpose of this function is to avoid re-running a
      complete T-factorization when a column is added. The function
      ONLY adds the column if it improves the rank of the matrix
      (nothing is done otherwise).  The function returns 1 if the
      column was added, and 0 otherwise.

      Note: T, U, P and the rank value supplied to the function must
      be correct one. This is not checked by the function (for reasons
      of efficiency). Also, the matrix operating on must have full
      column rank (this is assumed by the function, but not checked
      for reasons of efficiency).

      This function was written primarily intended for use by the LDPC
      codec class.
    */
    int T_fact_update_addcol(GF2mat &T, GF2mat &U, ivec &P, bvec newcol);

    // ----- Operators -----------

    //! Assignment operator
    void operator=(const GF2mat &X);

    //! Check if equal
    bool operator==(const GF2mat &X) const; 

    // ----- Friends ------

    //! Multiplication operator
    friend GF2mat operator*(const GF2mat &X, const GF2mat &Y);

    //! Multiplication operator with binary vector
    friend bvec operator*(const GF2mat &X, const bvec &y);

    /*! \brief Addition operator 

    Subtraction is not implemented because it is
      equivalent to addition. */
    friend GF2mat operator+(const GF2mat &X, const GF2mat &Y);

    //! Output stream operator (plain text)
    friend std::ostream &operator<<(std::ostream &os, const GF2mat &X);

    //! Write the matrix to file
    friend it_file &operator<<(it_file &f, const GF2mat &X);

    //! Read the matrix from file
    friend it_ifile &operator>>(it_ifile &f, GF2mat &X);

    //!Multiplication X*Y' where X and Y are GF(2) matrices
    friend GF2mat mult_trans(const GF2mat &X, const GF2mat &Y);

  private:
    int nwords;                  // number of bytes used
    int nrows, ncols;            // number of rows and columns of matrix
    Mat<unsigned short int> data;     // data structure
  };


  // ============ RELATED FUNCTIONS ============================

  /*!  
    /relates GF2mat 
    /brief Write GF(2) matrix to file.
  */
  it_file &operator<<(it_file &f, const GF2mat &X);

  /*!  
    /relates GF2mat 
    /brief Read GF(2) matrix from file.
  */
  it_ifile &operator>>(it_ifile &f, GF2mat &X);

  /*!
    \relates GF2mat
    \brief GF(2) matrix multiplication 
  */
  GF2mat operator*(const GF2mat &X, const GF2mat &Y);

  /*!
    \relates GF2mat
    \brief GF(2) matrix multiplication with "regular" binary vector
  */
  bvec operator*(const GF2mat &X, const bvec &y);

  /*!
    \relates GF2mat
    \brief GF(2) matrix addition 
  */
  GF2mat operator+(const GF2mat &X, const GF2mat &Y);

  /*!
    \relates GF2mat
    \brief Output stream (plain text) operator for dense GF(2) matrices
  */
  std::ostream &operator<<(std::ostream &os, const GF2mat &X);

  /*!
    \relates GF2mat
    \brief GF(2) Identity matrix
  */
  GF2mat gf2dense_eye(int m);

  /*! \relates GF2mat
    \brief Multiplication X*Y' where X and Y are GF(2) matrices
  */
  GF2mat mult_trans(const GF2mat &X, const GF2mat &Y);

  // ================= INLINE IMPLEMENTATIONS =======================

  inline void GF2mat::addto_element(int i, int j, bin s) 
  {
    it_assert0(i>=0 && i<nrows,"GF2mat::addto_element()");
    it_assert0(j>=0 && j<ncols,"GF2mat::addto_element()");
    if (s==1) { 
      int col = j>>lImax;
      char bit = j&mImax;
      data(i,col) ^= (1<<bit); 
    }
  }

  inline bin GF2mat::get(int i, int j) const 
  {
    it_assert0(i>=0 && i<nrows,"GF2mat::get_element()");
    it_assert0(j>=0 && j<ncols,"GF2mat::get_element()");
    int col = j>>lImax;
    char bit = j&mImax;
    return ((data(i,col) >> bit) & 1);    // NB data must be unsigned for this to be well defined
  }

  inline void GF2mat::set(int i, int j, bin s) 
  {
    it_assert0(i>=0 && i<nrows,"GF2mat::set_element()");
    it_assert0(j>=0 && j<ncols,"GF2mat::set_element()");
    it_assert0(s==0 || s==1,"GF2mat::set_element()");
    int col = j>>lImax;
    char bit = j&mImax;
    //    int oldvalue = (data(i,col) >> bit) & 1;
    if (s==1) {  // set bit to one
      data(i,col) |=  (1<<bit);
    }  else { // set bit to zero
      data(i,col) &= ((~0) ^ (1<<bit));
    }
  }

} // namespace itpp

#endif // #ifndef GF2MAT_H
