/*!
 * \file
 * \brief Matrix Class Definitions
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

#ifndef MAT_H
#define MAT_H

#ifndef _MSC_VER
#  include <itpp/config.h>
#else
#  include <itpp/config_msvc.h>
#endif

#include <itpp/base/itassert.h>
#include <itpp/base/math/misc.h>
#include <itpp/base/factory.h>


namespace itpp {

  // Declaration of Vec
  template<class Num_T> class Vec;
  // Declaration of Mat
  template<class Num_T> class Mat;
  // Declaration of bin
  class bin;

  //! Horizontal concatenation of two matrices
  template<class Num_T> const Mat<Num_T> concat_horizontal(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
  //! Vertical concatenation of two matrices
  template<class Num_T> const Mat<Num_T> concat_vertical(const Mat<Num_T> &m1, const Mat<Num_T> &m2);

  //! Addition of two matricies
  template<class Num_T> const Mat<Num_T> operator+(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
  //! Addition of a matrix and a scalar
  template<class Num_T> const Mat<Num_T> operator+(const Mat<Num_T> &m, Num_T t);
  //! Addition of a scalar and a matrix
  template<class Num_T> const Mat<Num_T> operator+(Num_T t, const Mat<Num_T> &m);

  //! Subtraction of two matrices
  template<class Num_T> const Mat<Num_T> operator-(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
  //! Subtraction of matrix and scalar
  template<class Num_T> const Mat<Num_T> operator-(const Mat<Num_T> &m, Num_T t);
  //! Subtraction of scalar and matrix
  template<class Num_T> const Mat<Num_T> operator-(Num_T t, const Mat<Num_T> &m);
  //! Negation of matrix
  template<class Num_T> const Mat<Num_T> operator-(const Mat<Num_T> &m);

  //! Multiplication of two matricies
  template<class Num_T> const Mat<Num_T> operator*(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
  //! Multiplication of matrix and vector
  template<class Num_T> const Vec<Num_T> operator*(const Mat<Num_T> &m, const Vec<Num_T> &v);
  //! Multiplication of vector and matrix (only works for a matrix that is a row vector)
  template<class Num_T> const Mat<Num_T> operator*(const Vec<Num_T> &v, const Mat<Num_T> &m);
  //! Multiplication of matrix and scalar
  template<class Num_T> const Mat<Num_T> operator*(const Mat<Num_T> &m, Num_T t);
  //! Multiplication of scalar and matrix
  template<class Num_T> const Mat<Num_T> operator*(Num_T t, const Mat<Num_T> &m);

  //! Element wise multiplication of two matrices. Same functionality as Matlab/Octave expression A .* B
  template<class Num_T> const Mat<Num_T> elem_mult(const Mat<Num_T> &A, const Mat<Num_T> &B);
  //! Element wise multiplication of two matrices, storing the result in matrix \c out (which is re-sized if necessary)
  template<class Num_T> void elem_mult_out(const Mat<Num_T> &A, const Mat<Num_T> &B, Mat<Num_T> &out);
  //! Element wise multiplication of three matrices, storing the result in matrix \c out (which is re-sized if necessary)
  template<class Num_T> void elem_mult_out(const Mat<Num_T> &A, const Mat<Num_T> &B, const Mat<Num_T> &C, Mat<Num_T> &out);
  //! Element wise multiplication of four matrices, storing the result in matrix \c out (which is re-sized if necessary)
  template<class Num_T> void elem_mult_out(const Mat<Num_T> &A, const Mat<Num_T> &B, const Mat<Num_T> &C, const Mat<Num_T> &D, Mat<Num_T> &out);
  //! In-place element wise multiplication of two matrices. Fast version of B = elem_mult(A,B)
  template<class Num_T> void elem_mult_inplace(const Mat<Num_T> &A, Mat<Num_T> &B);
  //! Element wise multiplication of two matrices, followed by summation of the resultant elements. Fast version of sumsum(elem_mult(A,B))
  template<class Num_T> Num_T elem_mult_sum(const Mat<Num_T> &A, const Mat<Num_T> &B);

  //! Division of matrix and scalar
  template<class Num_T> const Mat<Num_T> operator/(const Mat<Num_T> &m, Num_T t);

  //! Element wise division of two matrices. Same functionality as Matlab/Octave expression A ./ B
  template<class Num_T> const Mat<Num_T> elem_div(const Mat<Num_T> &A, const Mat<Num_T> &B);
  //! Element wise division of two matrices, storing the result in matrix \c out (which is re-sized if necessary)
  template<class Num_T> void elem_div_out(const Mat<Num_T> &A, const Mat<Num_T> &B, Mat<Num_T> &out);
  //! Element wise division of two matrices, followed by summation of the resultant elements. Fast version of sumsum(elem_div(A,B))
  template<class Num_T> Num_T elem_div_sum(const Mat<Num_T> &A, const Mat<Num_T> &B);

  // -------------------------------------------------------------------------------------
  // Declaration of Mat
  // -------------------------------------------------------------------------------------

  /*!
    \ingroup arr_vec_mat
    \brief Matrix Class (Templated)
    \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson

    Matrices can be of arbitrarily types, but conversions and functions are
    prepared for \c bin, \c short, \c int, \c double, and \c complex<double>
    vectors and these are predefined as: \c bmat, \c smat, \c imat, \c mat,
    and \c cmat. \c double and \c complex<double> are usually \c double and
    \c complex<double> respectively. However, this can be changed when
    compiling the it++ (see installation notes for more details). (Note: for
    binary matrices, an alternative to the bmat class is \c GF2mat and
    \c GF2mat_dense, which offer a more memory efficient representation and
    additional functions for linear algebra.)

    Examples:

    Matrix Constructors:
    When constructing a matrix without dimensions (memory) use
    \code mat temp; \endcode
    For construction of a matrix of a given size use
    \code mat temp(rows, cols); \endcode
    It is also possible to assign the constructed matrix the value and dimension
    of another matrix by
    \code vec temp(inmatrix); \endcode
    If you have explicit values you would like to assign to the matrix it is
    possible to do this using strings as:
    \code
    mat a("0 0.7;5 9.3"); // that is a = [0, 0.7; 5, 9.3]
    mat a="0 0.7;5 9.3";  // the constructor are called implicitly
    \endcode
    It is also possible to change dimension by
    \code temp.set_size(new_rows, new_cols, false); \endcode
    where \c false is used to indicate that the old values in \c temp
    is not copied. If you like to preserve the values use \c true.

    There are a number of methods to access parts of a matrix. Examples are
    \code
    a(5,3);     // Element number (5,3)
    a(5,9,3,5);  // Sub-matrix from rows 5, 6, 7, 8, 9 the columns 3, 4, and 5
    a.get_row(10);  // Row 10
    a.get_col(10); // Column 10
    \endcode

    It is also possible to modify parts of a vector as e.g. in
    \code
    a.set_row(5, invector);    // Set row 5 to \c invector
    a.set_col(3, invector); // Set column 3 to \c invector
    a.copy_col(1, 5); // Copy column 5 to column 1
    a.swap_cols(1, 5); // Swap the contents of columns 1 and 5
    \endcode

    It is of course also possible to perform the common linear algebra
    methods such as addition, subtraction, and matrix multiplication. Observe
    though, that vectors are assumed to be column-vectors in operations with
    matrices.

    Most elementary functions such as sin(), cosh(), log(), abs(), ..., are
    available as operations on the individual elements of the matrices. Please
    see the individual functions for more details.

    By default, the Mat elements are created using the default constructor for
    the element type. This can be changed by specifying a suitable Factory in
    the Mat constructor call; see Detailed Description for Factory.
  */
  template<class Num_T>
  class Mat {
  public:
    //! The type of the matrix values
    typedef Num_T value_type;

    //! Default constructor. An element factory \c f can be specified
    explicit Mat(const Factory &f = DEFAULT_FACTORY);
    //! Create a matrix of size (inrow, incol). An element factory \c f can be specified
    Mat(int inrow, int incol, const Factory &f = DEFAULT_FACTORY);
    //! Copy constructor
    Mat(const Mat<Num_T> &m);
    //! Constructor, similar to the copy constructor, but also takes an element factory \c f as argument
    Mat(const Mat<Num_T> &m, const Factory &f);
    //! Create a copy of the vector \c invector treated as a column vector. An element factory \c f can be specified
    Mat(const Vec<Num_T> &invector, const Factory &f = DEFAULT_FACTORY);
    //! Set matrix equal to values in string. An element factory \c f can be specified
    Mat(const std::string &str, const Factory &f = DEFAULT_FACTORY);
    //! Set matrix equal to values in string. An element factory \c f can be specified
    Mat(const char *str, const Factory &f = DEFAULT_FACTORY);
    /*!
      \brief Constructor taking a C-array as input. Copies all data. An element factory \c f can be specified

      By default the matrix is stored as a RowMajor matrix (i.e. listing elements in sequence
      beginning with the first column).
    */
    Mat(const Num_T *c_array, int rows, int cols, bool RowMajor = true, const Factory &f = DEFAULT_FACTORY);

    //! Destructor
    ~Mat();

    //! The number of columns
    int cols() const { return no_cols; }
    //! The number of rows
    int rows() const { return no_rows; }
    //! The number of elements
    int size() const { return datasize; }
    //! Set size of matrix. If copy = true then keep the data before resizing.
    void set_size(int inrow, int incol, bool copy = false);
    //! Set matrix equal to the all zero matrix
    void zeros();
    //! Set matrix equal to the all zero matrix
    void clear() { zeros(); }
    //! Set matrix equal to the all one matrix
    void ones();
    //! Set matrix equal to values in \c values
    bool set(const char *values);
    //! Set matrix equal to values in the string \c str
    bool set(const std::string &str);

    //! Get element (R,C) from matrix
    const Num_T &operator()(int R, int C) const;
    //! Get element (R,C) from matrix
    Num_T &operator()(int R, int C);
    //! Get element \c index using linear addressing (by rows)
    const Num_T &operator()(int index) const;
    //! Get element \c index using linear addressing (by rows)
    Num_T &operator()(int index);
    //! Get element (R,C) from matrix
    const Num_T &get(int R,int C) const;
    //! Set element (R,C) of matrix
    void set(int R,int C, const Num_T &v);

    /*!
      \brief Sub-matrix from row \c r1 to row \c r2 and columns \c c1 to \c c2.

      Value -1 indicates the last row and column, respectively.
    */
    const Mat<Num_T> operator()(int r1, int r2, int c1, int c2) const;
    /*!
      \brief Sub-matrix from row \c r1 to row \c r2 and columns \c c1 to \c c2.

      Value -1 indicates the last row and column, respectively.
    */
    const Mat<Num_T> get(int r1, int r2, int c1, int c2) const;

    //! Get row \c Index
    const Vec<Num_T> get_row(int Index) const ;
    //! Get rows \c r1 through \c r2
    const Mat<Num_T> get_rows(int r1, int r2) const;
    //! Get the rows specified by \c indexlist
    const Mat<Num_T> get_rows(const Vec<int> &indexlist) const;
    //! Get column \c Index
    const Vec<Num_T> get_col(int Index) const ;
    //! Get columns \c c1 through \c c2
    const Mat<Num_T> get_cols(int c1, int c2) const;
    //! Get the columns specified by \c indexlist
    const Mat<Num_T> get_cols(const Vec<int> &indexlist) const;
    //! Set row \c Index to \c invector
    void set_row(int Index, const Vec<Num_T> &invector);
    //! Set column \c Index to \c invector
    void set_col(int Index, const Vec<Num_T> &invector);
    //! Copy row \c from onto row \c to
    void copy_row(int to, int from);
    //! Copy column \c from onto column \c to
    void copy_col(int to, int from);
    //! Swap the rows \c r1 and \c r2
    void swap_rows(int r1, int r2);
    //! Swap the columns \c c1 and \c c2
    void swap_cols(int c1, int c2);

    //! Set submatrix defined by rows r1,r2 and columns c1,c2 to matrix m
    void set_submatrix(int r1, int r2, int c1, int c2, const Mat<Num_T> &m);
    //! Set submatrix defined by upper-left element (r,c) and the size of matrix m to m
    void set_submatrix(int r, int c, const Mat<Num_T> &m);
    //! Set all elements of submatrix defined by rows r1,r2 and columns c1,c2 to value t
    void set_submatrix(int r1, int r2, int c1, int c2, const Num_T t);

    //! Delete row number \c r
    void del_row(int r);
    //! Delete rows from \c r1 to \c r2
    void del_rows(int r1, int r2);
    //! Delete column number \c c
    void del_col(int c);
    //! Delete columns from \c c1 to \c c2
    void del_cols(int c1, int c2);
    //! Insert vector \c in at row number \c r, the matrix can be empty.
    void ins_row(int r, const Vec<Num_T> &in);
    //! Insert vector \c in at column number \c c, the matrix can be empty.
    void ins_col(int c, const Vec<Num_T> &in);
    //! Append vector \c to the bottom of the matrix, the matrix can be empty.
    void append_row(const Vec<Num_T> &in);
    //! Append vector \c to the right of the matrix, the matrix can be empty.
    void append_col(const Vec<Num_T> &in);

    //! Matrix transpose
    const Mat<Num_T> transpose() const;
    //! Matrix transpose
    const Mat<Num_T> T() const { return this->transpose(); }
    //! Hermitian matrix transpose (conjugate transpose)
    const Mat<Num_T> hermitian_transpose() const;
    //! Hermitian matrix transpose (conjugate transpose)
    const Mat<Num_T> H() const { return this->hermitian_transpose(); }

    //! Concatenate the matrices \c m1 and \c m2 horizontally
    friend const Mat<Num_T> concat_horizontal <>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
    //! Concatenate the matrices \c m1 and \c m2 vertically
    friend const Mat<Num_T> concat_vertical <>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);

    //! Set all elements of the matrix equal to \c t
    Mat<Num_T>& operator=(Num_T t);
    //! Set matrix equal to \c m
    Mat<Num_T>& operator=(const Mat<Num_T> &m);
    //! Set matrix equal to the vector \c v, assuming column vector
    Mat<Num_T>& operator=(const Vec<Num_T> &v);
    //! Set matrix equal to values in the string
    Mat<Num_T>& operator=(const char *values);

    //! Addition of matrices
    Mat<Num_T>& operator+=(const Mat<Num_T> &m);
    //! Addition of scalar to matrix
    Mat<Num_T>& operator+=(Num_T t);
    //! Addition of two matrices
    friend const  Mat<Num_T> operator+<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
    //! Addition of matrix and scalar
    friend const Mat<Num_T> operator+<>(const Mat<Num_T> &m, Num_T t);
    //! Addition of scalar and matrix
    friend const Mat<Num_T> operator+<>(Num_T t, const Mat<Num_T> &m);

    //! Subtraction of matrix
    Mat<Num_T>& operator-=(const Mat<Num_T> &m);
    //! Subtraction of scalar from matrix
    Mat<Num_T>& operator-=(Num_T t);
    //! Subtraction of \c m2 from \c m1
    friend const Mat<Num_T> operator-<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
    //! Subraction of scalar from matrix
    friend const Mat<Num_T> operator-<>(const Mat<Num_T> &m, Num_T t);
    //! Subtract matrix from scalar
    friend const Mat<Num_T> operator-<>(Num_T t, const Mat<Num_T> &m);
    //! Subraction of matrix
    friend const Mat<Num_T> operator-<>(const Mat<Num_T> &m);

    //! Matrix multiplication
    Mat<Num_T>& operator*=(const Mat<Num_T> &m);
    //! Multiplication by a scalar
    Mat<Num_T>& operator*=(Num_T t);
    //! Multiplication of two matrices
    friend const Mat<Num_T> operator*<>(const Mat<Num_T> &m1, const Mat<Num_T> &m2);
    //! Multiplication of matrix \c m and vector \c v (column vector)
    friend const Vec<Num_T> operator*<>(const Mat<Num_T> &m, const Vec<Num_T> &v);
    //! Multiplication of transposed vector \c v and matrix \c m
    friend const Mat<Num_T> operator*<>(const Vec<Num_T> &v, const Mat<Num_T> &m);
    //! Multiplication of matrix and scalar
    friend const Mat<Num_T> operator*<>(const Mat<Num_T> &m, Num_T t);
    //! Multiplication of scalar and matrix
    friend const Mat<Num_T> operator*<>(Num_T t, const Mat<Num_T> &m);

    //! Element wise multiplication of two matrices. Same functionality as Matlab expression A .* B
    friend const Mat<Num_T> elem_mult <>(const Mat<Num_T> &A, const Mat<Num_T> &B);
    //! Element wise multiplication of two matrices, storing the result in matrix \c out (which is re-sized if necessary)
    friend void elem_mult_out <>(const Mat<Num_T> &A, const Mat<Num_T> &B, Mat<Num_T> &out);
    //! Element wise multiplication of three matrices, storing the result in matrix \c out (which is re-sized if necessary)
    friend void elem_mult_out <>(const Mat<Num_T> &A, const Mat<Num_T> &B, const Mat<Num_T> &C, Mat<Num_T> &out);
    //! Element wise multiplication of four matrices, storing the result in matrix \c out (which is re-sized if necessary)
    friend void elem_mult_out <>(const Mat<Num_T> &A, const Mat<Num_T> &B, const Mat<Num_T> &C, const Mat<Num_T> &D, Mat<Num_T> &out);
    //! In-place element wise multiplication of two matrices. Fast version of B = elem_mult(A,B)
    friend void elem_mult_inplace <>(const Mat<Num_T> &A, Mat<Num_T> &B);
    //! Element wise multiplication of two matrices, followed by summation of the resultant elements. Fast version of sumsum(elem_mult(A,B))
    friend Num_T elem_mult_sum <>(const Mat<Num_T> &A, const Mat<Num_T> &B);

    //! Division by a scalar
    Mat<Num_T>& operator/=(Num_T t);
    //! Division of matrix with scalar
    friend const Mat<Num_T> operator/<>(const Mat<Num_T> &m, Num_T t);
    //! Elementwise division with the current matrix
    Mat<Num_T>& operator/=(const Mat<Num_T> &m);

    //! Element wise division of two matrices. Same functionality as Matlab expression A ./ B
    friend const Mat<Num_T> elem_div <>(const Mat<Num_T> &A, const Mat<Num_T> &B);
    //! Element wise division of two matrices, storing the result in matrix \c out (which is re-sized if necessary)
    friend void elem_div_out <>(const Mat<Num_T> &A, const Mat<Num_T> &B, Mat<Num_T> &out);
    //! Element wise division of two matrices, followed by summation of the resultant elements. Fast version of sumsum(elem_div(A,B))
    friend Num_T elem_div_sum <>(const Mat<Num_T> &A, const Mat<Num_T> &B);

    //! Compare two matrices. False if wrong sizes or different values
    bool operator==(const Mat<Num_T> &m) const;
    //! Compare two matrices. True if different
    bool operator!=(const Mat<Num_T> &m) const;

    //! Get element (R,C) from matrix without boundary check (Not recommended to use).
    Num_T &_elem(int R,int C) { return data[R+C*no_rows]; }
    //! Get element (R,C) from matrix without boundary check (Not recommended to use).
    const Num_T &_elem(int R,int C) const { return data[R+C*no_rows]; }
    //! Get element \c index using linear addressing (by rows) without boundary check (Not recommended to use).
    Num_T &_elem(int index) { return data[index]; }
    //! Get element \c index using linear addressing (by rows) without boundary check (Not recommended to use).
    const Num_T &_elem(int index) const { return data[index]; }

    //! Access of the internal data structure. Don't use. May be changed!
    Num_T *_data() { return data; }
    //! Access to the internal data structure. Don't use. May be changed!
    const Num_T *_data() const { return data; }
    //! Access to the internal data structure. Don't use. May be changed!
    int _datasize() const { return datasize; }

  protected:
    //! Allocate memory for the matrix
    void alloc(int rows, int cols);
    //! Free the memory space of the matrix
    void free();

    /*! Protected integer variables
     * @{ */
    int datasize, no_rows, no_cols;
    /*! @} */
    //! Protected data pointer
    Num_T *data;
    //! Element factory (set to DEFAULT_FACTORY to use Num_T default constructors only)
    const Factory &factory;
  };

  // -------------------------------------------------------------------------------------
  // Type definitions of mat, cmat, imat, smat, and bmat
  // -------------------------------------------------------------------------------------

  /*!
    \relates Mat
    \brief Default Matrix Type
  */
  typedef Mat<double> mat;

  /*!
    \relates Mat
    \brief Default Complex Matrix Type
  */
  typedef Mat<std::complex<double> > cmat;

  /*!
    \relates Mat
    \brief Integer matrix
  */
  typedef Mat<int> imat;

  /*!
    \relates Mat
    \brief short int matrix
  */
  typedef Mat<short int> smat;

  /*!
    \relates Mat
    \relates GF2mat
    \relates GF2mat_sparse
    \brief bin matrix
  */
  typedef Mat<bin> bmat;

} //namespace itpp


#include <itpp/base/vec.h>

namespace itpp {

  // ----------------------------------------------------------------------
  // Declaration of input and output streams for Mat
  // ----------------------------------------------------------------------

  /*!
    \relatesalso Mat
    \brief Output stream for matrices
  */
  template <class Num_T>
  std::ostream &operator<<(std::ostream &os, const Mat<Num_T> &m);

  /*!
    \relatesalso Mat
    \brief Input stream for matrices

    The input can be on the form "1 2 3; 4 5 6" or "[[1 2 3][4 5 6]]", i.e. with
    brackets or semicolons as row delimiters. The first form is compatible with
    the set method, while the second form is compatible with the ostream
    operator. The elements on a row can be separated by blank space or commas.
    Rows that are shorter than the longest row are padded with zero elements.
    "[]" means an empty matrix.
  */
  template <class Num_T>
  std::istream &operator>>(std::istream &is, Mat<Num_T> &m);

  // ----------------------------------------------------------------------
  // Implementation of templated Mat members and friends
  // ----------------------------------------------------------------------

  template<class Num_T> inline
  void Mat<Num_T>::alloc(int rows, int cols)
  {
//     if (datasize == rows * cols) { // Reuse the memory
//       no_rows = rows; no_cols = cols;
//       return;
//     }
//     free();  // Free memory (if any allocated)
    if ((rows > 0) && (cols > 0)) {
      datasize = rows * cols;
      no_rows = rows;
      no_cols = cols;
      create_elements(data, datasize, factory);
    }
    else {
      data = 0;
      datasize = 0;
      no_rows = 0;
      no_cols = 0;
    }
  }

  template<class Num_T> inline
  void Mat<Num_T>::free()
  {
    delete[] data;
    data = 0;
    datasize = 0;
    no_rows = 0;
    no_cols = 0;
  }


  template<class Num_T> inline
  Mat<Num_T>::Mat(const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f) {}

  template<class Num_T> inline
  Mat<Num_T>::Mat(int inrow, int incol, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
  {
    it_assert_debug((inrow >= 0) && (incol >= 0), "The rows and columns must be >= 0");
    alloc(inrow, incol);
  }

  template<class Num_T> inline
  Mat<Num_T>::Mat(const Mat<Num_T> &m) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(m.factory)
  {
    alloc(m.no_rows, m.no_cols);
    copy_vector(m.datasize, m.data, data);
  }

  template<class Num_T> inline
  Mat<Num_T>::Mat(const Mat<Num_T> &m, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
  {
    alloc(m.no_rows, m.no_cols);
    copy_vector(m.datasize, m.data, data);
  }

  template<class Num_T> inline
  Mat<Num_T>::Mat(const Vec<Num_T> &invector, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
  {
    int size = invector.size();
    alloc(size, 1);
    copy_vector(size, invector._data(), data);
  }

  template<class Num_T> inline
  Mat<Num_T>::Mat(const std::string &str, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
  {
    set(str);
  }

  template<class Num_T> inline
  Mat<Num_T>::Mat(const char *str, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
  {
    set(str);
  }

  template<class Num_T> inline
  Mat<Num_T>::Mat(const Num_T *c_array, int rows, int cols, bool RowMajor, const Factory &f) :
    datasize(0), no_rows(0), no_cols(0), data(0), factory(f)
  {
    alloc(rows, cols);

    if (!RowMajor)
      copy_vector(datasize, c_array, data);
    else { // Row Major
      for (int i=0; i<rows; i++) {
	for (int j=0; j<cols; j++)
	  data[i+j*no_rows] = c_array[i*no_cols+j];
      }
    }
  }

  template<class Num_T> inline
  Mat<Num_T>::~Mat()
  {
    free();
  }


  template<class Num_T> inline
  void Mat<Num_T>::set_size(int inrow, int incol, bool copy)
  {
    it_assert_debug((inrow >= 0) && (incol >= 0), "Mat<Num_T>::set_size: "
		    "The number of rows and columns must be >= 0");
    // check if we have to resize the current matrix
    if ((no_rows == inrow) && (no_cols == incol))
      return;
    // conditionally copy previous matrix content
    if (copy) {
      // backup access to the current matrix
      Num_T* tmp = data;
      int tmp_rows = no_rows;
      // check the boundaries of the copied data
      int min_r = (no_rows < inrow) ? no_rows : inrow;
      int min_c = (no_cols < incol) ? no_cols : incol;
      // allocate new memory
      alloc(inrow, incol);
      // copy the previous data into the allocated memory
      for (int i = 0; i < min_c; ++i) {
	copy_vector(min_r, &tmp[i*tmp_rows], &data[i*no_rows]);
      }
      // fill-in the rest of matrix with zeros
      for (int i = min_r; i < inrow; ++i)
	for (int j = 0; j < incol; ++j)
	  data[i+j*inrow] = Num_T(0);
      for (int j = min_c; j < incol; ++j)
	for (int i = 0; i < min_r; ++i)
	  data[i+j*inrow] = Num_T(0);
      // free the previous data memory
      delete[] tmp;
    }
    // if possible, reuse the allocated memory
    else if (datasize == inrow * incol) {
      no_rows = inrow;
      no_cols = incol;
    }
    // finally release old memory and allocate a new one
    else {
      free();
      alloc(inrow, incol);
    }
  }

  template<class Num_T> inline
  void Mat<Num_T>::zeros()
  {
    for(int i=0; i<datasize; i++)
      data[i] = Num_T(0);
  }

  template<class Num_T> inline
  void Mat<Num_T>::ones()
  {
    for(int i=0; i<datasize; i++)
      data[i] = Num_T(1);
  }

  template<class Num_T> inline
  const Num_T& Mat<Num_T>::operator()(int R,int C) const
  {
    it_assert_debug((R >= 0) && (R < no_rows) && (C >= 0) && (C < no_cols),
	       "Mat<Num_T>::operator(): index out of range");
    return data[R+C*no_rows];
  }

  template<class Num_T> inline
  Num_T& Mat<Num_T>::operator()(int R,int C)
  {
    it_assert_debug((R >= 0) && (R < no_rows) && (C >= 0) && (C < no_cols),
	       "Mat<Num_T>::operator(): index out of range");
    return data[R+C*no_rows];
  }

  template<class Num_T> inline
  Num_T& Mat<Num_T>::operator()(int index)
  {
    it_assert_debug((index < no_rows*no_cols) && (index >= 0),
	       "Mat<Num_T>::operator(): index out of range");
    return data[index];
  }

  template<class Num_T> inline
  const Num_T& Mat<Num_T>::operator()(int index) const
  {
    it_assert_debug((index < no_rows*no_cols) && (index >= 0),
	       "Mat<Num_T>::operator(): index out of range");
    return data[index];
  }

  template<class Num_T> inline
  const Num_T& Mat<Num_T>::get(int R,int C) const
  {
    it_assert_debug((R >= 0) && (R < no_rows) && (C >= 0) && (C < no_cols),
	       "Mat<Num_T>::get(): index out of range");
    return data[R+C*no_rows];
  }

  template<class Num_T> inline
  void Mat<Num_T>::set(int R, int C, const Num_T &v)
  {
    it_assert_debug((R >= 0) && (R < no_rows) && (C >= 0) && (C < no_cols),
	       "Mat<Num_T>::set(): index out of range");
    data[R+C*no_rows] = v;
  }


  template<class Num_T>
  bool Mat<Num_T>::set(const char *values)
  {
    std::istringstream buffer(values);
    int rows = 0, maxrows = 10, cols = 0, nocols = 0, maxcols = 10;
    bool comma = true;

    free();
    alloc(maxrows, maxcols);
    zeros();

    while (buffer.peek() != EOF) {
      if (++rows > maxrows) {
	maxrows <<= 1;
	set_size(maxrows, maxcols, true);
      }

      cols = 0;
      while ((buffer.peek() != ';') && (buffer.peek() != EOF)) {
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

	default:
	  if (++cols > nocols) {
	    nocols = cols;
	    if (cols > maxcols) {
	      maxcols <<= 1;
	      set_size(maxrows, maxcols, true);
	    }
	  }
	  buffer >> this->operator()(rows-1, cols-1);
	  it_assert(!buffer.fail(), "Mat<Num_T>::set(): Stream operation failed (buffer >> data)");
	  comma = false;
	}
      }

      if (buffer.peek() == ';') {
	// remove row separator...
	buffer.seekg(1, std::ios_base::cur);
	// ... but also check if it was at the end of buffer
	while (buffer.peek() == ' ' || buffer.peek() == '\t') {
	  buffer.seekg(1, std::ios_base::cur);
	}
	comma = false;
	// it_assert(buffer.peek() != EOF, "Mat<double>::set(): Improper data string");
      }

      it_assert(!comma || (nocols == 0), "Mat<Num_T>::set(): Improper matrix string");

    }
    set_size(rows, nocols, true);

    return true;
  }

  //! Specialization of \c set() method for int
  template<>
  bool Mat<int>::set(const char *values);
  //! Specialization of \c set() method for short int
  template<>
  bool Mat<short int>::set(const char *values);

  template<class Num_T>
  bool Mat<Num_T>::set(const std::string &str)
  {
    return set(str.c_str());
  }

  template<class Num_T> inline
  const Mat<Num_T> Mat<Num_T>::operator()(int r1, int r2, int c1, int c2) const
  {
    if (r1 == -1) r1 = no_rows-1;
    if (r2 == -1) r2 = no_rows-1;
    if (c1 == -1) c1 = no_cols-1;
    if (c2 == -1) c2 = no_cols-1;

    it_assert_debug(r1>=0 && r2>=0 && r1<no_rows && r2<no_rows &&
	       c1>=0 && c2>=0 && c1<no_cols && c2<no_cols, "operator()(r1,r2,c1,c2)");

    it_assert_debug(r2>=r1 && c2>=c1, "Mat<Num_T>::op(): r2>=r1 or c2>=c1 not fulfilled");

    Mat<Num_T> s(r2-r1+1, c2-c1+1);

    for (int i=0;i<s.no_cols;i++)
      copy_vector(s.no_rows, data+r1+(c1+i)*no_rows, s.data+i*s.no_rows);

    return s;
  }

  template<class Num_T> inline
  const Mat<Num_T> Mat<Num_T>::get(int r1, int r2, int c1, int c2) const
  {
    return (*this)(r1, r2, c1, c2);
  }

  template<class Num_T> inline
  const Vec<Num_T> Mat<Num_T>::get_row(int Index) const
  {
    it_assert_debug(Index>=0 && Index<no_rows, "Mat<Num_T>::get_row: index out of range");
    Vec<Num_T> a(no_cols);

    copy_vector(no_cols, data+Index, no_rows, a._data(), 1);
    return a;
  }

  template<class Num_T>
  const Mat<Num_T> Mat<Num_T>::get_rows(int r1, int r2) const
  {
    it_assert_debug(r1>=0 && r2<no_rows && r1<=r2, "Mat<Num_T>::get_rows: index out of range");
    Mat<Num_T> m(r2-r1+1, no_cols);

    for (int i=0; i<m.rows(); i++)
      copy_vector(no_cols, data+i+r1, no_rows, m.data+i, m.no_rows);

    return m;
  }

  template<class Num_T>
  const Mat<Num_T> Mat<Num_T>::get_rows(const Vec<int> &indexlist) const
  {
    Mat<Num_T> m(indexlist.size(),no_cols);

    for (int i=0;i<indexlist.size();i++) {
      it_assert_debug(indexlist(i)>=0 && indexlist(i)<no_rows, "Mat<Num_T>::get_rows: index out of range");
      copy_vector(no_cols, data+indexlist(i), no_rows, m.data+i, m.no_rows);
    }

    return m;
  }

  template<class Num_T> inline
  const Vec<Num_T> Mat<Num_T>::get_col(int Index) const
  {
    it_assert_debug(Index>=0 && Index<no_cols, "Mat<Num_T>::get_col: index out of range");
    Vec<Num_T> a(no_rows);

    copy_vector(no_rows, data+Index*no_rows, a._data());

    return a;
  }

  template<class Num_T> inline
  const Mat<Num_T> Mat<Num_T>::get_cols(int c1, int c2) const
  {
    it_assert_debug(c1>=0 && c2<no_cols && c1<=c2, "Mat<Num_T>::get_cols: index out of range");
    Mat<Num_T> m(no_rows, c2-c1+1);

    for (int i=0; i<m.cols(); i++)
      copy_vector(no_rows, data+(i+c1)*no_rows, m.data+i*m.no_rows);

    return m;
  }

  template<class Num_T> inline
  const Mat<Num_T> Mat<Num_T>::get_cols(const Vec<int> &indexlist) const
  {
    Mat<Num_T> m(no_rows,indexlist.size());

    for (int i=0; i<indexlist.size(); i++) {
      it_assert_debug(indexlist(i)>=0 && indexlist(i)<no_cols, "Mat<Num_T>::get_cols: index out of range");
      copy_vector(no_rows, data+indexlist(i)*no_rows, m.data+i*m.no_rows);
    }

    return m;
  }

  template<class Num_T> inline
  void Mat<Num_T>::set_row(int Index, const Vec<Num_T> &v)
  {
    it_assert_debug(Index>=0 && Index<no_rows, "Mat<Num_T>::set_row: index out of range");
    it_assert_debug(v.length() == no_cols, "Mat<Num_T>::set_row: lengths doesn't match");

    copy_vector(v.size(), v._data(), 1, data+Index, no_rows);
  }

  template<class Num_T> inline
  void Mat<Num_T>::set_col(int Index, const Vec<Num_T> &v)
  {
    it_assert_debug(Index>=0 && Index<no_cols, "Mat<Num_T>::set_col: index out of range");
    it_assert_debug(v.length() == no_rows, "Mat<Num_T>::set_col: lengths doesn't match");

    copy_vector(v.size(), v._data(), data+Index*no_rows);
  }

  template<class Num_T> inline
  void Mat<Num_T>::copy_row(int to, int from)
  {
    it_assert_debug(to>=0 && from>=0 && to<no_rows && from<no_rows,
	       "Mat<Num_T>::copy_row: index out of range");

    if (from == to)
      return;

    copy_vector(no_cols, data+from, no_rows, data+to, no_rows);
  }

  template<class Num_T> inline
  void Mat<Num_T>::copy_col(int to, int from)
  {
    it_assert_debug(to>=0 && from>=0 && to<no_cols && from<no_cols,
	       "Mat<Num_T>::copy_col: index out of range");

    if (from == to)
      return;

    copy_vector(no_rows, data+from*no_rows, data+to*no_rows);
  }

  template<class Num_T> inline
  void Mat<Num_T>::swap_rows(int r1, int r2)
  {
    it_assert_debug(r1>=0 && r2>=0 && r1<no_rows && r2<no_rows,
	       "Mat<Num_T>::swap_rows: index out of range");

    if (r1 == r2)
      return;

    swap_vector(no_cols, data+r1, no_rows, data+r2, no_rows);
  }

  template<class Num_T> inline
  void Mat<Num_T>::swap_cols(int c1, int c2)
  {
    it_assert_debug(c1>=0 && c2>=0 && c1<no_cols && c2<no_cols,
	       "Mat<Num_T>::swap_cols: index out of range");

    if (c1 == c2)
      return;

    swap_vector(no_rows, data+c1*no_rows, data+c2*no_rows);
  }

  template<class Num_T> inline
  void Mat<Num_T>::set_submatrix(int r1, int r2, int c1, int c2, const Mat<Num_T> &m)
  {

    if (r1 == -1) r1 = no_rows-1;
    if (r2 == -1) r2 = no_rows-1;
    if (c1 == -1) c1 = no_cols-1;
    if (c2 == -1) c2 = no_cols-1;

    it_assert_debug(r1>=0 && r2>=0 && r1<no_rows && r2<no_rows &&
	       c1>=0 && c2>=0 && c1<no_cols && c2<no_cols, "Mat<Num_T>::set_submatrix(): index out of range");

    it_assert_debug(r2>=r1 && c2>=c1, "Mat<Num_T>::set_submatrix: r2<r1 or c2<c1");
    it_assert_debug(m.no_rows == r2-r1+1 && m.no_cols == c2-c1+1, "Mat<Num_T>::set_submatrix(): sizes don't match");

    for (int i=0; i<m.no_cols; i++)
      copy_vector(m.no_rows, m.data+i*m.no_rows, data+(c1+i)*no_rows+r1);
  }



  template<class Num_T> inline
  void Mat<Num_T>::set_submatrix(int r, int c, const Mat<Num_T> &m)
  {

    it_assert_debug(r>=0 && r+m.no_rows<=no_rows &&
	       c>=0 && c+m.no_cols<=no_cols, "Mat<Num_T>::set_submatrix(): index out of range");

    for (int i=0; i<m.no_cols; i++)
      copy_vector(m.no_rows, m.data+i*m.no_rows, data+(c+i)*no_rows+r);
  }



  template<class Num_T> inline
  void Mat<Num_T>::set_submatrix(int r1, int r2, int c1, int c2, const Num_T t)
  {

    if (r1 == -1) r1 = no_rows-1;
    if (r2 == -1) r2 = no_rows-1;
    if (c1 == -1) c1 = no_cols-1;
    if (c2 == -1) c2 = no_cols-1;

    it_assert_debug(r1>=0 && r2>=0 && r1<no_rows && r2<no_rows &&
	       c1>=0 && c2>=0 && c1<no_cols && c2<no_cols, "Mat<Num_T>::set_submatrix(): index out of range");

    it_assert_debug(r2>=r1 && c2>=c1, "Mat<Num_T>::set_submatrix: r2<r1 or c2<c1");

    int i, j, pos, rows = r2-r1+1;

    for (i=c1; i<=c2; i++) {
      pos = i*no_rows+r1;
      for (j=0; j<rows; j++) {
	data[pos++] = t;
      }
    }
  }

  template<class Num_T> inline
  void Mat<Num_T>::del_row(int r)
  {
    it_assert_debug(r>=0 && r<no_rows, "Mat<Num_T>::del_row(): index out of range");
    Mat<Num_T> Temp(*this);
    set_size(no_rows-1, no_cols, false);
    for (int i=0 ; i < r ; i++) {
      copy_vector(no_cols, &Temp.data[i], no_rows+1, &data[i], no_rows);
    }
    for (int i=r ; i < no_rows ; i++) {
      copy_vector(no_cols, &Temp.data[i+1], no_rows+1, &data[i], no_rows);
    }

  }

  template<class Num_T> inline
  void Mat<Num_T>::del_rows(int r1, int r2)
  {
    it_assert_debug((r1 >= 0) && (r2 < no_rows) && (r1 <= r2),
		    "Mat<Num_T>::del_rows(): indices out of range");

    Mat<Num_T> Temp(*this);
    int no_del_rows = r2-r1+1;
    set_size(no_rows-no_del_rows, no_cols, false);
    for (int i = 0; i < r1 ; ++i) {
      copy_vector(no_cols, &Temp.data[i], Temp.no_rows, &data[i], no_rows);
    }
    for (int i = r2+1; i < Temp.no_rows; ++i) {
      copy_vector(no_cols, &Temp.data[i], Temp.no_rows, &data[i-no_del_rows],
		  no_rows);
    }
  }

  template<class Num_T> inline
  void Mat<Num_T>::del_col(int c)
  {
    it_assert_debug(c>=0 && c<no_cols, "Mat<Num_T>::del_col(): index out of range");
    Mat<Num_T> Temp(*this);

    set_size(no_rows, no_cols-1, false);
    copy_vector(c*no_rows, Temp.data, data);
    copy_vector((no_cols - c)*no_rows, &Temp.data[(c+1)*no_rows], &data[c*no_rows]);
  }

  template<class Num_T> inline
  void Mat<Num_T>::del_cols(int c1, int c2)
  {
    it_assert_debug(c1>=0 && c2<no_cols && c1<=c2, "Mat<Num_T>::del_cols(): index out of range");
    Mat<Num_T> Temp(*this);
    int n_deleted_cols = c2-c1+1;
    set_size(no_rows, no_cols-n_deleted_cols, false);
    copy_vector(c1*no_rows, Temp.data, data);
    copy_vector((no_cols-c1)*no_rows, &Temp.data[(c2+1)*no_rows], &data[c1*no_rows]);
  }

  template<class Num_T> inline
  void Mat<Num_T>::ins_row(int r, const Vec<Num_T> &in)
  {
    it_assert_debug(r>=0 && r<=no_rows, "Mat<Num_T>::ins_row(): index out of range");
    it_assert_debug((in.size() == no_cols) || no_rows==0, "Mat<Num_T>::ins_row(): vector size does not match");

    if (no_cols==0) {
      no_cols = in.size();
    }

    Mat<Num_T> Temp(*this);
    set_size(no_rows+1, no_cols, false);

    for (int i=0 ; i < r ; i++) {
      copy_vector(no_cols, &Temp.data[i], no_rows-1, &data[i], no_rows);
    }
    copy_vector(no_cols, in._data(), 1, &data[r], no_rows);
    for (int i=r+1 ; i < no_rows ; i++) {
      copy_vector(no_cols, &Temp.data[i-1], no_rows-1, &data[i], no_rows);
    }
  }

  template<class Num_T> inline
  void Mat<Num_T>::ins_col(int c, const Vec<Num_T> &in)
  {
    it_assert_debug(c>=0 && c<=no_cols, "Mat<Num_T>::ins_col(): index out of range");
    it_assert_debug(in.size() == no_rows || no_cols==0, "Mat<Num_T>::ins_col(): vector size does not match");

    if (no_rows==0) {
      no_rows = in.size();
    }

    Mat<Num_T> Temp(*this);
    set_size(no_rows, no_cols+1, false);

    copy_vector(c*no_rows, Temp.data, data);
    copy_vector(no_rows, in._data(), &data[c*no_rows]);
    copy_vector((no_cols-c-1)*no_rows, &Temp.data[c*no_rows], &data[(c+1)*no_rows]);
  }

  template<class Num_T> inline
  void Mat<Num_T>::append_row(const Vec<Num_T> &in)
  {
    ins_row(no_rows, in);
  }

  template<class Num_T> inline
  void Mat<Num_T>::append_col(const Vec<Num_T> &in)
  {
    ins_col(no_cols, in);
  }

  template<class Num_T>
  const Mat<Num_T> Mat<Num_T>::transpose() const
  {
    Mat<Num_T> temp(no_cols, no_rows);
    for (int i = 0; i < no_rows; ++i) {
      copy_vector(no_cols, &data[i], no_rows, &temp.data[i * no_cols], 1);
    }
    return temp;
  }

  //! Specialization of \c hermitian_transpose() method for std::complex<double>
  template<>
  const cmat Mat<std::complex<double> >::hermitian_transpose() const;

  template<class Num_T>
  const Mat<Num_T> Mat<Num_T>::hermitian_transpose() const
  {
    Mat<Num_T> temp(no_cols, no_rows);
    for (int i = 0; i < no_rows; ++i) {
      copy_vector(no_cols, &data[i], no_rows, &temp.data[i * no_cols], 1);
    }
    return temp;
  }

  template<class Num_T>
  const Mat<Num_T> concat_horizontal(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
  {
    it_assert_debug(m1.no_rows == m2.no_rows,
		    "Mat<Num_T>::concat_horizontal(): wrong sizes");
    int no_rows = m1.no_rows;
    Mat<Num_T> temp(no_rows, m1.no_cols + m2.no_cols);
    for (int i = 0; i < m1.no_cols; ++i) {
      copy_vector(no_rows, &m1.data[i * no_rows], &temp.data[i * no_rows]);
    }
    for (int i = 0; i < m2.no_cols; ++i) {
      copy_vector(no_rows, &m2.data[i * no_rows], &temp.data[(m1.no_cols + i)
							     * no_rows]);
    }
    return temp;
  }

  template<class Num_T>
  const Mat<Num_T> concat_vertical(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
  {
    it_assert_debug(m1.no_cols == m2.no_cols,
		    "Mat<Num_T>::concat_vertical; wrong sizes");
    int no_cols = m1.no_cols;
    Mat<Num_T> temp(m1.no_rows + m2.no_rows, no_cols);
    for (int i = 0; i < no_cols; ++i) {
      copy_vector(m1.no_rows, &m1.data[i * m1.no_rows],
		  &temp.data[i * temp.no_rows]);
      copy_vector(m2.no_rows, &m2.data[i * m2.no_rows],
		  &temp.data[i * temp.no_rows + m1.no_rows]);
    }
    return temp;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator=(Num_T t)
  {
    for (int i=0; i<datasize; i++)
      data[i] = t;
    return *this;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator=(const Mat<Num_T> &m)
  {
    if (this != &m) {
      set_size(m.no_rows,m.no_cols, false);
      if (m.datasize != 0)
	copy_vector(m.datasize, m.data, data);
    }
    return *this;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator=(const Vec<Num_T> &v)
  {
    it_assert_debug((no_rows == 1 && no_cols == v.size())
		    || (no_cols == 1 && no_rows == v.size()),
		    "Mat<Num_T>::operator=(): Wrong size of the argument");
    set_size(v.size(), 1, false);
    copy_vector(v.size(), v._data(), data);
    return *this;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator=(const char *values)
  {
    set(values);
    return *this;
  }

  //-------------------- Templated friend functions --------------------------

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator+=(const Mat<Num_T> &m)
  {
    if (datasize == 0)
      operator=(m);
    else {
      int i, j, m_pos=0, pos=0;
      it_assert_debug(m.no_rows==no_rows && m.no_cols==no_cols,"Mat<Num_T>::operator+=: wrong sizes");
      for (i=0; i<no_cols; i++) {
	for (j=0; j<no_rows; j++)
	  data[pos+j] += m.data[m_pos+j];
	pos += no_rows;
	m_pos += m.no_rows;
      }
    }
    return *this;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator+=(Num_T t)
  {
    for (int i=0; i<datasize; i++)
      data[i] += t;
    return *this;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator+(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
  {
    Mat<Num_T> r(m1.no_rows, m1.no_cols);
    int i, j, m1_pos=0, m2_pos=0, r_pos=0;

    it_assert_debug(m1.no_rows==m2.no_rows && m1.no_cols == m2.no_cols, "Mat<Num_T>::operator+: wrong sizes");

    for (i=0; i<r.no_cols; i++) {
      for (j=0; j<r.no_rows; j++)
	r.data[r_pos+j] = m1.data[m1_pos+j] + m2.data[m2_pos+j];
      // next column
      m1_pos += m1.no_rows;
      m2_pos += m2.no_rows;
      r_pos += r.no_rows;
    }

    return r;
  }


  template<class Num_T> inline
  const Mat<Num_T> operator+(const Mat<Num_T> &m, Num_T t)
  {
    Mat<Num_T> r(m.no_rows, m.no_cols);

    for (int i=0; i<r.datasize; i++)
      r.data[i] = m.data[i] + t;

    return r;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator+(Num_T t, const Mat<Num_T> &m)
  {
    Mat<Num_T> r(m.no_rows, m.no_cols);

    for (int i=0; i<r.datasize; i++)
      r.data[i] = t + m.data[i];

    return r;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator-=(const Mat<Num_T> &m)
  {
    int i,j, m_pos=0, pos=0;

    if (datasize == 0) {
      set_size(m.no_rows, m.no_cols, false);
      for (i=0; i<no_cols; i++) {
	for (j=0; j<no_rows; j++)
	  data[pos+j] = -m.data[m_pos+j];
	// next column
	m_pos += m.no_rows;
	pos += no_rows;
      }
    } else {
      it_assert_debug(m.no_rows==no_rows && m.no_cols==no_cols,"Mat<Num_T>::operator-=: wrong sizes");
      for (i=0; i<no_cols; i++) {
	for (j=0; j<no_rows; j++)
	  data[pos+j] -= m.data[m_pos+j];
	// next column
	m_pos += m.no_rows;
	pos += no_rows;
      }
    }
    return *this;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator-(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
  {
    Mat<Num_T> r(m1.no_rows, m1.no_cols);
    int i, j, m1_pos=0, m2_pos=0, r_pos=0;

    it_assert_debug(m1.no_rows==m2.no_rows && m1.no_cols == m2.no_cols, "Mat<Num_T>::operator-: wrong sizes");

    for (i=0; i<r.no_cols; i++) {
      for (j=0; j<r.no_rows; j++)
	r.data[r_pos+j] = m1.data[m1_pos+j] - m2.data[m2_pos+j];
      // next column
      m1_pos += m1.no_rows;
      m2_pos += m2.no_rows;
      r_pos += r.no_rows;
    }

    return r;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator-=(Num_T t)
  {
    for (int i=0; i<datasize; i++)
      data[i] -= t;
    return *this;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator-(const Mat<Num_T> &m, Num_T t)
  {
    Mat<Num_T> r(m.no_rows, m.no_cols);
    int i, j, m_pos=0, r_pos=0;

    for (i=0; i<r.no_cols; i++) {
      for (j=0; j<r.no_rows; j++)
	r.data[r_pos+j] = m.data[m_pos+j] - t;
      // next column
      m_pos += m.no_rows;
      r_pos += r.no_rows;
    }

    return r;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator-(Num_T t, const Mat<Num_T> &m)
  {
    Mat<Num_T> r(m.no_rows, m.no_cols);
    int i, j, m_pos=0, r_pos=0;

    for (i=0; i<r.no_cols; i++) {
      for (j=0; j<r.no_rows; j++)
	r.data[r_pos+j] = t - m.data[m_pos+j];
      // next column
      m_pos += m.no_rows;
      r_pos += r.no_rows;
    }

    return r;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator-(const Mat<Num_T> &m)
  {
    Mat<Num_T> r(m.no_rows, m.no_cols);
    int i, j, m_pos=0, r_pos=0;

    for (i=0; i<r.no_cols; i++) {
      for (j=0; j<r.no_rows; j++)
	r.data[r_pos+j] = -m.data[m_pos+j];
      // next column
      m_pos += m.no_rows;
      r_pos += r.no_rows;
    }

    return r;
  }

#if defined(HAVE_CBLAS)
  template<> mat& Mat<double>::operator*=(const mat &m);
  template<> cmat& Mat<std::complex<double> >::operator*=(const cmat &m);
#endif

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator*=(const Mat<Num_T> &m)
  {
    it_assert_debug(no_cols == m.no_rows,"Mat<Num_T>::operator*=: wrong sizes");
    Mat<Num_T> r(no_rows, m.no_cols);

    Num_T tmp;

    int i,j,k, r_pos=0, pos=0, m_pos=0;

    for (i=0; i<r.no_cols; i++) {
      for (j=0; j<r.no_rows; j++) {
	tmp = Num_T(0);
	pos = 0;
	for (k=0; k<no_cols; k++) {
	  tmp += data[pos+j] * m.data[m_pos+k];
	  pos += no_rows;
	}
	r.data[r_pos+j] = tmp;
      }
      r_pos += r.no_rows;
      m_pos += m.no_rows;
    }
    operator=(r); // time consuming
    return *this;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator*=(Num_T t)
  {
    for (int i=0; i<datasize; i++)
      data[i] *= t;
    return *this;
  }

#if defined(HAVE_CBLAS)
  template<> const mat operator*(const mat &m1, const mat &m2);
  template<> const cmat operator*(const cmat &m1, const cmat &m2);
#endif


  template<class Num_T>
  const Mat<Num_T> operator*(const Mat<Num_T> &m1, const Mat<Num_T> &m2)
  {
    it_assert_debug(m1.no_cols == m2.no_rows,"Mat<Num_T>::operator*: wrong sizes");
    Mat<Num_T> r(m1.no_rows, m2.no_cols);

    Num_T tmp;
    int i, j, k;
    Num_T *tr=r.data, *t1, *t2=m2.data;

    for (i=0; i<r.no_cols; i++) {
      for (j=0; j<r.no_rows; j++) {
	tmp = Num_T(0); t1 = m1.data+j;
	for (k=m1.no_cols; k>0; k--) {
	  tmp += *(t1) * *(t2++);
	  t1 += m1.no_rows;
	}
	*(tr++) = tmp; t2 -= m2.no_rows;
      }
      t2 += m2.no_rows;
    }

    return r;
  }

#if defined(HAVE_CBLAS)
  template<> const vec operator*(const mat &m, const vec &v);
  template<> const cvec operator*(const cmat &m, const cvec &v);
#endif

  template<class Num_T>
  const Vec<Num_T> operator*(const Mat<Num_T> &m, const Vec<Num_T> &v)
  {
    it_assert_debug(m.no_cols == v.size(),"Mat<Num_T>::operator*: wrong sizes");
    Vec<Num_T> r(m.no_rows);
    int i, k, m_pos;

    for (i=0; i<m.no_rows; i++) {
      r(i) = Num_T(0);
      m_pos = 0;
      for (k=0; k<m.no_cols; k++) {
	r(i) += m.data[m_pos+i] * v(k);
	m_pos += m.no_rows;
      }
    }

    return r;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator*(const Vec<Num_T> &v, const Mat<Num_T> &m)
  {
    it_assert_debug(m.no_rows == 1,"Mat<Num_T>::operator*: wrong sizes");
    Mat<Num_T> r(v.size(), m.no_cols);

    for (int i = 0; i < v.size(); ++i)
      for (int j = 0; j < m.no_cols; ++j)
	r(i,j) = v(i) * m.data[j];

    return r;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator*(const Mat<Num_T> &m, Num_T t)
  {
    Mat<Num_T> r(m.no_rows, m.no_cols);

    for (int i=0; i<r.datasize; i++)
      r.data[i] = m.data[i] * t;

    return r;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator*(Num_T t, const Mat<Num_T> &m)
  {
    Mat<Num_T> r(m.no_rows, m.no_cols);

    for (int i=0; i<r.datasize; i++)
      r.data[i] = m.data[i] * t;

    return r;
  }

  template<class Num_T> inline
  const Mat<Num_T> elem_mult(const Mat<Num_T> &A, const Mat<Num_T> &B)
  {
    Mat<Num_T> out;
    elem_mult_out(A,B,out);
    return out;
  }

  template<class Num_T> inline
  void elem_mult_out(const Mat<Num_T> &A, const Mat<Num_T> &B, Mat<Num_T> &out)
  {
    it_assert_debug( (A.no_rows==B.no_rows) && (A.no_cols==B.no_cols), "Mat<Num_T>::elem_mult_out: wrong sizes");

    if( (out.no_rows != A.no_rows) || (out.no_cols != A.no_cols) )
      out.set_size(A.no_rows, A.no_cols);

    for(int i=0; i<out.datasize; i++)
      out.data[i] = A.data[i] * B.data[i];
  }

  template<class Num_T> inline
  void elem_mult_out(const Mat<Num_T> &A, const Mat<Num_T> &B, const Mat<Num_T> &C, Mat<Num_T> &out)
  {
    it_assert_debug( (A.no_rows==B.no_rows==C.no_rows) \
                && (A.no_cols==B.no_cols==C.no_cols), \
                "Mat<Num_T>::elem_mult_out: wrong sizes" );

    if( (out.no_rows != A.no_rows) || (out.no_cols != A.no_cols) )
      out.set_size(A.no_rows, A.no_cols);

    for(int i=0; i<out.datasize; i++)
      out.data[i] = A.data[i] * B.data[i] * C.data[i];
  }

  template<class Num_T> inline
  void elem_mult_out(const Mat<Num_T> &A, const Mat<Num_T> &B, const Mat<Num_T> &C, const Mat<Num_T> &D, Mat<Num_T> &out)
  {
    it_assert_debug( (A.no_rows==B.no_rows==C.no_rows==D.no_rows) \
                && (A.no_cols==B.no_cols==C.no_cols==D.no_cols), \
                "Mat<Num_T>::elem_mult_out: wrong sizes" );
    if( (out.no_rows != A.no_rows) || (out.no_cols != A.no_cols) )
      out.set_size(A.no_rows, A.no_cols);

    for(int i=0; i<out.datasize; i++)
      out.data[i] = A.data[i] * B.data[i] * C.data[i] * D.data[i];
  }

  template<class Num_T> inline
  void elem_mult_inplace(const Mat<Num_T> &A, Mat<Num_T> &B)
  {
    it_assert_debug( (A.no_rows==B.no_rows) && (A.no_cols==B.no_cols), \
                "Mat<Num_T>::elem_mult_inplace: wrong sizes" );

    for(int i=0; i<B.datasize; i++)
      B.data[i] *= A.data[i];
  }

  template<class Num_T> inline
  Num_T elem_mult_sum(const Mat<Num_T> &A, const Mat<Num_T> &B)
  {
    it_assert_debug( (A.no_rows==B.no_rows) && (A.no_cols==B.no_cols), "Mat<Num_T>::elem_mult_sum: wrong sizes" );

    Num_T acc = 0;

    for(int i=0; i<A.datasize; i++)
      acc += A.data[i] * B.data[i];

    return acc;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator/=(Num_T t)
  {
    for (int i=0; i<datasize; i++)
      data[i] /= t;
    return *this;
  }

  template<class Num_T> inline
  const Mat<Num_T> operator/(const Mat<Num_T> &m, Num_T t)
  {
    Mat<Num_T> r(m.no_rows, m.no_cols);

    for (int i=0; i<r.datasize; i++)
      r.data[i] = m.data[i] / t;

    return r;
  }

  template<class Num_T> inline
  Mat<Num_T>& Mat<Num_T>::operator/=(const Mat<Num_T> &m)
  {
    it_assert_debug(m.no_rows==no_rows && m.no_cols==no_cols, "Mat<Num_T>::operator/=: wrong sizes");

    for (int i=0; i<datasize; i++)
      data[i] /= m.data[i];
    return *this;
  }

  template<class Num_T> inline
  const Mat<Num_T> elem_div(const Mat<Num_T> &A, const Mat<Num_T> &B)
  {
    Mat<Num_T> out;
    elem_div_out(A,B,out);
    return out;
  }

  template<class Num_T> inline
  void elem_div_out(const Mat<Num_T> &A, const Mat<Num_T> &B, Mat<Num_T> &out)
  {
    it_assert_debug( (A.no_rows==B.no_rows) && (A.no_cols==B.no_cols), "Mat<Num_T>::elem_div_out: wrong sizes");

    if( (out.no_rows != A.no_rows) || (out.no_cols != A.no_cols) )
      out.set_size(A.no_rows, A.no_cols);

    for(int i=0; i<out.datasize; i++)
      out.data[i] = A.data[i] / B.data[i];
  }

  template<class Num_T> inline
  Num_T elem_div_sum(const Mat<Num_T> &A, const Mat<Num_T> &B)
  {
    it_assert_debug( (A.no_rows==B.no_rows) && (A.no_cols==B.no_cols), "Mat<Num_T>::elem_div_sum: wrong sizes" );

    Num_T acc = 0;

    for(int i=0; i<A.datasize; i++)
      acc += A.data[i] / B.data[i];

    return acc;
  }

  template<class Num_T>
  bool Mat<Num_T>::operator==(const Mat<Num_T> &m) const
  {
    if (no_rows!=m.no_rows || no_cols != m.no_cols) return false;
    for (int i=0;i<datasize;i++) {
      if (data[i]!=m.data[i]) return false;
    }
    return true;
  }

  template<class Num_T>
  bool Mat<Num_T>::operator!=(const Mat<Num_T> &m) const
  {
    if (no_rows != m.no_rows || no_cols != m.no_cols) return true;
    for (int i=0;i<datasize;i++) {
      if (data[i]!=m.data[i]) return true;
    }
    return false;
  }

  template <class Num_T>
  std::ostream &operator<<(std::ostream &os, const Mat<Num_T> &m)
  {
    int i;

    switch (m.rows()) {
    case 0 :
      os << "[]";
      break;
    case 1 :
      os << '[' << m.get_row(0) << ']';
      break;
    default:
      os << '[' << m.get_row(0) << std::endl;
      for (i=1; i<m.rows()-1; i++)
	      os << ' ' << m.get_row(i) << std::endl;
      os << ' ' << m.get_row(m.rows()-1) << ']';
    }

    return os;
  }

  template <class Num_T>
  std::istream &operator>>(std::istream &is, Mat<Num_T> &m)
  {
    std::ostringstream buffer;
    bool started = false;
    bool finished = false;
    bool brackets = false;
    bool within_double_brackets = false;
    char c;

    while (!finished) {
      if (is.eof()) {
        finished = true;
      } else {
        c = is.get();

        if (is.eof() || (c == '\n')) {
          if (brackets) {
            // Right bracket missing
            is.setstate(std::ios_base::failbit);
            finished = true;
          } else if (!((c == '\n') && !started)) {
            finished = true;
          }
        } else if ((c == ' ') || (c == '\t')) {
          if (started) {
            buffer << ' ';
          }
        } else if (c == '[') {
          if ((started && !brackets) || within_double_brackets) {
            // Unexpected left bracket
            is.setstate(std::ios_base::failbit);
            finished = true;
          } else if (!started) {
            started = true;
            brackets = true;
          } else {
            within_double_brackets = true;
          }
        } else if (c == ']') {
          if (!started || !brackets) {
            // Unexpected right bracket
            is.setstate(std::ios_base::failbit);
            finished = true;
          } else if (within_double_brackets) {
            within_double_brackets = false;
            buffer << ';';
          } else {
            finished = true;
          }
          while (!is.eof() && (((c = is.peek()) == ' ') || (c == '\t'))) {
            is.get();
          }
          if (!is.eof() && (c == '\n')) {
            is.get();
          }
        } else {
          started = true;
          buffer << c;
        }
      }
    }

    if (!started) {
      m.set_size(0, false);
    } else {
      m.set(buffer.str());
    }

    return is;
  }

#ifndef _MSC_VER

  //---------------------------------------------------------------------
  // Instantiations
  //---------------------------------------------------------------------

  //------- class instantiations --------

  //! Template instantiation of Mat<double>
  extern template class Mat<double>;
  //! Template instantiation of Mat<std::complex<double> >
  extern template class Mat<std::complex<double> >;
  //! Template instantiation of Mat<int>
  extern template class Mat<int>;
  //! Template instantiation of Mat<short int>
  extern template class Mat<short int>;
  //! Template instantiation of Mat<bin>
  extern template class Mat<bin>;


  //-------------------- Operator instantiations --------------------

  //-------- Addition operators ---------------

  //! Template instantiation of operator+
  extern template const mat operator+(const mat &m1, const mat &m2);
  //! Template instantiation of operator+
  extern template const cmat operator+(const cmat &m1, const cmat &m2);
  //! Template instantiation of operator+
  extern template const imat operator+(const imat &m1, const imat &m2);
  //! Template instantiation of operator+
  extern template const smat operator+(const smat &m1, const smat &m2);
  //! Template instantiation of operator+
  extern template const bmat operator+(const bmat &m1, const bmat &m2);

  //! Template instantiation of operator+
  extern template const mat operator+(const mat &m, double t);
  //! Template instantiation of operator+
  extern template const cmat operator+(const cmat &m, std::complex<double> t);
  //! Template instantiation of operator+
  extern template const imat operator+(const imat &m, int t);
  //! Template instantiation of operator+
  extern template const smat operator+(const smat &m, short t);
  //! Template instantiation of operator+
  extern template const bmat operator+(const bmat &m, bin t);

  //! Template instantiation of operator+
  extern template const mat operator+(double t, const mat &m);
  //! Template instantiation of operator+
  extern template const cmat operator+(std::complex<double> t, const cmat &m);
  //! Template instantiation of operator+
  extern template const imat operator+(int t, const imat &m);
  //! Template instantiation of operator+
  extern template const smat operator+(short t, const smat &m);
  //! Template instantiation of operator+
  extern template const bmat operator+(bin t, const bmat &m);

  //-------- Subraction operators ---------------

  //! Template instantiation of operator-
  extern template const mat operator-(const mat &m1, const mat &m2);
  //! Template instantiation of operator-
  extern template const cmat operator-(const cmat &m1, const cmat &m2);
  //! Template instantiation of operator-
  extern template const imat operator-(const imat &m1, const imat &m2);
  //! Template instantiation of operator-
  extern template const smat operator-(const smat &m1, const smat &m2);
  //! Template instantiation of operator-
  extern template const bmat operator-(const bmat &m1, const bmat &m2);

  //! Template instantiation of operator-
  extern template const mat operator-(const mat &m, double t);
  //! Template instantiation of operator-
  extern template const cmat operator-(const cmat &m, std::complex<double> t);
  //! Template instantiation of operator-
  extern template const imat operator-(const imat &m, int t);
  //! Template instantiation of operator-
  extern template const smat operator-(const smat &m, short t);
  //! Template instantiation of operator-
  extern template const bmat operator-(const bmat &m, bin t);

  //! Template instantiation of operator-
  extern template const mat operator-(double t, const mat &m);
  //! Template instantiation of operator-
  extern template const cmat operator-(std::complex<double> t, const cmat &m);
  //! Template instantiation of operator-
  extern template const imat operator-(int t, const imat &m);
  //! Template instantiation of operator-
  extern template const smat operator-(short t, const smat &m);
  //! Template instantiation of operator-
  extern template const bmat operator-(bin t, const bmat &m);

  //--------- Unary minus ---------------

  //! Template instantiation of operator-
  extern template const mat operator-(const mat &m);
  //! Template instantiation of operator-
  extern template const cmat operator-(const cmat &m);
  //! Template instantiation of operator-
  extern template const imat operator-(const imat &m);
  //! Template instantiation of operator-
  extern template const smat operator-(const smat &m);
  //! Template instantiation of operator-
  extern template const bmat operator-(const bmat &m);

  //-------- Multiplication operators ---------------

#if !defined(HAVE_CBLAS)
  //! Template instantiation of operator*
  extern template const mat operator*(const mat &m1, const mat &m2);
  //! Template instantiation of operator*
  extern template const cmat operator*(const cmat &m1, const cmat &m2);
#endif
  //! Template instantiation of operator*
  extern template const imat operator*(const imat &m1, const imat &m2);
  //! Template instantiation of operator*
  extern template const smat operator*(const smat &m1, const smat &m2);
  //! Template instantiation of operator*
  extern template const bmat operator*(const bmat &m1, const bmat &m2);

#if !defined(HAVE_CBLAS)
  //! Template instantiation of operator*
  extern template const vec operator*(const mat &m, const vec &v);
  //! Template instantiation of operator*
  extern template const cvec operator*(const cmat &m, const cvec &v);
#endif
  //! Template instantiation of operator*
  extern template const ivec operator*(const imat &m, const ivec &v);
  //! Template instantiation of operator*
  extern template const svec operator*(const smat &m, const svec &v);
  //! Template instantiation of operator*
  extern template const bvec operator*(const bmat &m, const bvec &v);

  //! Template instantiation of operator*
  extern template const mat operator*(const vec &v, const mat &m);
  //! Template instantiation of operator*
  extern template const cmat operator*(const cvec &v, const cmat &m);
  //! Template instantiation of operator*
  extern template const imat operator*(const ivec &v, const imat &m);
  //! Template instantiation of operator*
  extern template const smat operator*(const svec &v, const smat &m);
  //! Template instantiation of operator*
  extern template const bmat operator*(const bvec &v, const bmat &m);

  //! Template instantiation of operator*
  extern template const mat operator*(const mat &m, double t);
  //! Template instantiation of operator*
  extern template const cmat operator*(const cmat &m, std::complex<double> t);
  //! Template instantiation of operator*
  extern template const imat operator*(const imat &m, int t);
  //! Template instantiation of operator*
  extern template const smat operator*(const smat &m, short t);
  //! Template instantiation of operator*
  extern template const bmat operator*(const bmat &m, bin t);

  //! Template instantiation of operator*
  extern template const mat operator*(double t, const mat &m);
  //! Template instantiation of operator*
  extern template const cmat operator*(std::complex<double> t, const cmat &m);
  //! Template instantiation of operator*
  extern template const imat operator*(int t, const imat &m);
  //! Template instantiation of operator*
  extern template const smat operator*(short t, const smat &m);
  //! Template instantiation of operator*
  extern template const bmat operator*(bin t, const bmat &m);

  // ------------ Elementwise multiplication -----------

  //! Template instantiation of elem_mult
  extern template const mat elem_mult(const mat &A, const mat &B);
  //! Template instantiation of elem_mult
  extern template const cmat elem_mult(const cmat &A, const cmat &B);
  //! Template instantiation of elem_mult
  extern template const imat elem_mult(const imat &A, const imat &B);
  // Extern Template Const instantiation of elem_mult
  //extern template const llmat elem_mult(const llmat &A, const llmat &B);
  //! Template instantiation of elem_mult
  extern template const smat elem_mult(const smat &A, const smat &B);
  //! Template instantiation of elem_mult
  extern template const bmat elem_mult(const bmat &A, const bmat &B);

  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const mat &A, const mat &B, mat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const cmat &A, const cmat &B, cmat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const imat &A, const imat &B, imat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const smat &A, const smat &B, smat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const bmat &A, const bmat &B, bmat &out);

  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const mat &A, const mat &B, const mat &C, mat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const cmat &A, const cmat &B, const cmat &C, cmat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const imat &A, const imat &B, const imat &C, imat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const smat &A, const smat &B, const smat &C, smat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const bmat &A, const bmat &B, const bmat &C, bmat &out);

  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const mat &A, const mat &B, const mat &C, const mat &D, mat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const cmat &A, const cmat &B, const cmat &C, const cmat &D, cmat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const imat &A, const imat &B, const imat &C, const imat &D, imat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const smat &A, const smat &B, const smat &C, const smat &D, smat &out);
  //! Template instantiation of elem_mult_out
  extern template void elem_mult_out(const bmat &A, const bmat &B, const bmat &C, const bmat &D, bmat &out);

  //! Template instantiation of elem_mult_inplace
  extern template void elem_mult_inplace(const mat &A, mat &B);
  //! Template instantiation of elem_mult_inplace
  extern template void elem_mult_inplace(const cmat &A, cmat &B);
  //! Template instantiation of elem_mult_inplace
  extern template void elem_mult_inplace(const imat &A, imat &B);
  //! Template instantiation of elem_mult_inplace
  extern template void elem_mult_inplace(const smat &A, smat &B);
  //! Template instantiation of elem_mult_inplace
  extern template void elem_mult_inplace(const bmat &A, bmat &B);

  //! Template instantiation of elem_mult_sum
  extern template double elem_mult_sum(const mat &A, const mat &B);
  //! Template instantiation of elem_mult_sum
  extern template std::complex<double> elem_mult_sum(const cmat &A, const cmat &B);
  //! Template instantiation of elem_mult_sum
  extern template int elem_mult_sum(const imat &A, const imat &B);
  //! Template instantiation of elem_mult_sum
  extern template short elem_mult_sum(const smat &A, const smat &B);
  //! Template instantiation of elem_mult_sum
  extern template bin elem_mult_sum(const bmat &A, const bmat &B);

  // ------------ Division operator -----------

  //! Template instantiation of operator/
  extern template const mat operator/(const mat &m, double t);
  //! Template instantiation of operator/
  extern template const cmat operator/(const cmat &m, std::complex<double> t);
  //! Template instantiation of operator/
  extern template const imat operator/(const imat &m, int t);
  //! Template instantiation of operator/
  extern template const smat operator/(const smat &m, short t);
  //! Template instantiation of operator/
  extern template const bmat operator/(const bmat &m, bin t);

  // ------------ Elementwise division -----------

  //! Template instantiation of elem_div
  extern template const mat elem_div(const mat &A, const mat &B);
  //! Template instantiation of elem_div
  extern template const cmat elem_div(const cmat &A, const cmat &B);
  //! Template instantiation of elem_div
  extern template const imat elem_div(const imat &A, const imat &B);
  //! Template instantiation of elem_div
  extern template const smat elem_div(const smat &A, const smat &B);
  //! Template instantiation of elem_div
  extern template const bmat elem_div(const bmat &A, const bmat &B);

  //! Template instantiation of elem_div_out
  extern template void elem_div_out(const mat &A, const mat &B, mat &out);
  //! Template instantiation of elem_div_out
  extern template void elem_div_out(const cmat &A, const cmat &B, cmat &out);
  //! Template instantiation of elem_div_out
  extern template void elem_div_out(const imat &A, const imat &B, imat &out);
  //! Template instantiation of elem_div_out
  extern template void elem_div_out(const smat &A, const smat &B, smat &out);
  //! Template instantiation of elem_div_out
  extern template void elem_div_out(const bmat &A, const bmat &B, bmat &out);

  //! Template instantiation of elem_div_sum
  extern template double elem_div_sum(const mat &A, const mat &B);
  //! Template instantiation of elem_div_sum
  extern template std::complex<double> elem_div_sum(const cmat &A, const cmat &B);
  //! Template instantiation of elem_div_sum
  extern template int elem_div_sum(const imat &A, const imat &B);
  //! Template instantiation of elem_div_sum
  extern template short elem_div_sum(const smat &A, const smat &B);
  //! Template instantiation of elem_div_sum
  extern template bin elem_div_sum(const bmat &A, const bmat &B);

  // ------------- Concatenations -----------------

  //! Template instantiation of concat_horizontal
  extern template const mat concat_horizontal(const mat &m1, const mat &m2);
  //! Template instantiation of concat_horizontal
  extern template const cmat concat_horizontal(const cmat &m1, const cmat &m2);
  //! Template instantiation of concat_horizontal
  extern template const imat concat_horizontal(const imat &m1, const imat &m2);
  //! Template instantiation of concat_horizontal
  extern template const smat concat_horizontal(const smat &m1, const smat &m2);
  //! Template instantiation of concat_horizontal
  extern template const bmat concat_horizontal(const bmat &m1, const bmat &m2);

  //! Template instantiation of concat_vertical
  extern template const mat concat_vertical(const mat &m1, const mat &m2);
  //! Template instantiation of concat_vertical
  extern template const cmat concat_vertical(const cmat &m1, const cmat &m2);
  //! Template instantiation of concat_vertical
  extern template const imat concat_vertical(const imat &m1, const imat &m2);
  //! Template instantiation of concat_vertical
  extern template const smat concat_vertical(const smat &m1, const smat &m2);
  //! Template instantiation of concat_vertical
  extern template const bmat concat_vertical(const bmat &m1, const bmat &m2);

  //----------- Output stream --------------

  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream &os, const mat  &m);
  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream &os, const cmat &m);
  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream &os, const imat  &m);
  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream &os, const smat  &m);
  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream &os, const bmat  &m);

  //----------- Input stream -------------

  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream &is, mat  &m);
  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream &is, cmat &m);
  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream &is, imat  &m);
  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream &is, smat  &m);
  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream &is, bmat  &m);

#endif // #ifndef _MSC_VER

} // namespace itpp

#endif // #ifndef MAT_H
