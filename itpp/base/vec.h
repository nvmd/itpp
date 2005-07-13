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
  \brief Templated Vector Class Definitions
  \author Tony Ottosson and Tobias Ringström

  $Revision$

  $Date$ 
*/

#ifndef __vec_h
#define __vec_h

#include <string> 
#include <iostream>
#include <complex>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <sstream>
#include "itpp/itconfig.h"
#include "itpp/base/binary.h"
#include "itpp/base/itassert.h"
#include "itpp/base/scalfunc.h"
#include "itpp/base/factory.h"
#include "itpp/base/copy_vector.h"



//using std::cout;
//using std::endl;
//using std::string;
//using std::ostream;
//using std::istream;
//using std::istringstream;
//using std::getline;
//using std::complex;

namespace itpp {

  // Declaration of Vec
  template<class Num_T> class Vec;
  // Declaration of Mat
  template<class Num_T> class Mat;
  // Declaration of bin
  class bin;

  //-----------------------------------------------------------------------------------
  // Declaration of Vec Friends 
  //-----------------------------------------------------------------------------------

  //! Addition of two vectors
  template<class Num_T> Vec<Num_T> operator+(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Addition of a vector and a scalar
  template<class Num_T> Vec<Num_T> operator+(const Vec<Num_T> &v, Num_T t);
  //! Addition of a scalar and a vector
  template<class Num_T> Vec<Num_T> operator+(Num_T t, const Vec<Num_T> &v);

  //! Subtraction of a vector from a vector
  template<class Num_T> Vec<Num_T> operator-(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Subtraction of a scalar from a vector
  template<class Num_T> Vec<Num_T> operator-(const Vec<Num_T> &v, Num_T t);
  //! Subtraction of vector from scalar. Results in a vector
  template<class Num_T> Vec<Num_T> operator-(Num_T t, const Vec<Num_T> &v);
  //! Negation of vector
  template<class Num_T> Vec<Num_T> operator-(const Vec<Num_T> &v);

  //! Inner (dot) product of two vectors v1 and v2
  template<class Num_T> Num_T dot(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Inner (dot) product of two vectors v1 and v2
  template<class Num_T> Num_T operator*(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  { return dot(v1, v2); }
  //! Outer product of two vectors v1 and v2
  template<class Num_T> Mat<Num_T> outer_product(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Multiplication of a vector and a scalar
  template<class Num_T> Vec<Num_T> operator*(const Vec<Num_T> &v, Num_T t);
  //! Multiplication of a scalar and a vector. Results in a vector
  template<class Num_T> Vec<Num_T> operator*(Num_T t, const Vec<Num_T> &v);
  //! Elementwise multiplication of the two vectors
  template<class Num_T> Vec<Num_T> elem_mult(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Elementwise multiplication of the three vectors
  template<class Num_T> Vec<Num_T> elem_mult(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3);
  //! Elementwise multiplication of the four vectors
  template<class Num_T> Vec<Num_T> elem_mult(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3, const Vec<Num_T> &v4);

  //! Division of all elements in \c v with \c t
  template<class Num_T> Vec<Num_T> operator/(const Vec<Num_T> &v, Num_T t);
  //! Division of \c t with all elements in \c v
  template<class Num_T> Vec<Num_T> operator/(const Num_T t, const Vec<Num_T> &v);
  //! Elementwise division
  template<class Num_T> Vec<Num_T> elem_div(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
  //! Elementwise division of scalar \c t and vector \c v
  template<class Num_T> Vec<Num_T> elem_div(const Num_T t, const Vec<Num_T> &v);

  //! Append element \c a to the end of the vector \c v
  template<class Num_T> Vec<Num_T> concat(const Vec<Num_T> &v, const Num_T a);
  //! Concat element \c a to the beginning of the vector \c v
  template<class Num_T> Vec<Num_T> concat(const Num_T a, const Vec<Num_T> &v);
  //! Concat vectors \c v1 and \c v2
  template<class Num_T> Vec<Num_T> concat(const Vec<Num_T> &v1,const Vec<Num_T> &v2);
  //! Concat vectors \c v1, \c v2 and \c v3
  template<class Num_T> Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3);
  //! Concat vectors \c v1, \c v2, \c v3 and \c v4
  template<class Num_T> Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3, const Vec<Num_T> &v4);
  //! Concat vectors \c v1, \c v2 \c v3, \c v4 and \c v5
  template<class Num_T> Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3,
					  const Vec<Num_T> &v4, const Vec<Num_T> &v5);

  //-----------------------------------------------------------------------------------
  // Declaration of Vec
  //-----------------------------------------------------------------------------------

  /*!
    \brief Templated vectors
    \author Tony Ottosson and Tobias Ringstrom

    Vectors can be of arbitrarily types, but conversions and functions are
    prepared for \c bin, \c short, \c int, \c double, and \c complex<double>
    vectors and these are predefined as: \c bvec, \c svec, \c ivec, \c vec,
    and \c cvec. \c double and \c complex<double> are \c double and
    \c complex<double> respectively.

    Examples:

    Vector Constructors:
    When constructing a vector without a length (memory) use
    \code vec temp; \endcode
    For construction of a vector of a given length use
    \code vec temp(length); \endcode
    It is also possible to assign the constructed vector the value and size
    of another vector by
    \code vec temp(invector); \endcode
    If you have explicit values you would like to assign to the vector it is
    possible to do this using strings as:
    \code
    vec a("0 0.7 5 9.3"); // that is a = [0, 0.7, 5, 9.3]
    vec a="0 0.7 5 9.3";  // the constructor are called implicitly
    ivec b="0:5";  // that is b = [0, 1, 2, 3, 4, 5]
    vec c="3:2.5:13";  // that is c = [3, 5.5, 8, 10.5, 13]
    \endcode
    It is also possible to change length by
    \code temp.set_size(new_length, false); \endcode
    where \c false is used to indicate that the old values in \c temp
    is not copied. If you like to preserve the values use \c true.

    There are a number of methods to access parts of a vector. Examples are
    \code
    a(5);     // Element number 5
    a(5,9);  // Elements 5, 6, 7, 8, and 9
    a.left(10);  // The 10 most left elements (the first)
    a.right(10); // The 10 most right elements (the last)
    a.mid(5, 7); // 7 elements starting from element 5
    \endcode

    It is also possible to modify parts of a vector as e.g. in
    \code
    a.del(5);    // deletes element number 5
    a.ins(3.4, 9); // inserts the element 3.4 at position 9
    a.replace_mid(12, b); // replaces elements from 12 with the vector b
    \endcode

    It is of course also possible to perform the common linear algebra
    methods such as addition, subtraction, and scalar product (*). Observe
    though, that vectors are assumed to be column-vectors in operations with
    matrices.

    Most elementary functions such as sin(), cosh(), log(), abs(), ..., are
    available as operations on the individual elements of the vectors. Please
    see the individual functions for more details.

    By default, the Vec elements are created using the default constructor for
    the element type. This can be changed by specifying a suitable Factory in
    the Vec constructor call; see Detailed Description for Factory.
  */
  template<class Num_T>
    class Vec {
  public:
    //! Default constructor. An element factory \c f can be specified
    explicit Vec(const Factory &f = DEFAULT_FACTORY) : factory(f) { init(); }
    //! Constructor. An element factory \c f can be specified
    explicit Vec(int size, const Factory &f = DEFAULT_FACTORY) : factory(f) { it_assert1(size>=0,"Negative size in Vec::Vec(int)"); init(); alloc(size); }
    //! Copy constructor
    Vec(const Vec<Num_T> &v);
    //! Constructor, similar to the copy constructor, but also takes an element factory \c f as argument
    Vec(const Vec<Num_T> &v, const Factory &f);
    //! Constructor. An element factory \c f can be specified
    Vec(const char *values, const Factory &f = DEFAULT_FACTORY) : factory(f) { init(); set(values); }
    //! Constructor. An element factory \c f can be specified
    Vec(const std::string &values, const Factory &f = DEFAULT_FACTORY) : factory(f) { init(); set(values); }
    //! Constructor taking a C-array as input. Copies all data. An element factory \c f can be specified
    Vec(Num_T *c_array, int size, const Factory &f = DEFAULT_FACTORY) : factory(f) { init(); alloc(size); copy_vector(size, c_array, data); }

    //! Destructor
    ~Vec() { free(); }

    //! The size of the vector
    int length() const { return datasize; }
    //! The size of the vector
    int size() const { return datasize; }

    //! Set length of vector. if copy = true then keeping the old values
    void set_length(int size, bool copy=false) { set_size(size,copy); }
    //! Set length of vector. if copy = true then keeping the old values
    void set_size(int size, bool copy=false);
    //! Set the vector to the all zero vector
    void zeros() { for (int i=0; i<datasize; i++) {data[i]=Num_T(0);} }
    //! Set the vector to the all zero vector
    void clear() { zeros(); }
    //! Set the vector to the all one vector
    void ones() { for (int i=0; i<datasize; i++) {data[i]=Num_T(1);} }
    //! Set the vector equal to the values in the \c str string
    bool set(const char *str);
    //! Set the vector equal to the values in the \c str string
    bool set(const std::string &str);

    //! C-style index operator. First element is 0
    const Num_T &operator[](int i) const { it_assert0(i>=0&&i<datasize, "operator[]"); return data[i]; }
    //! Index operator. First element is 0
    const Num_T &operator()(int i) const { it_assert0(i>=0&&i<datasize, "operator()"); return data[i]; }
    //! C-style index operator. First element is 0
    Num_T &operator[](int i) { it_assert0(i>=0&&i<datasize, "operator[]"); return data[i]; }
    //! Index operator. First element is 0
    Num_T &operator()(int i) { it_assert0(i>=0&&i<datasize, "operator()"); return data[i]; }
    //! Sub-vector with elements from \c i1 to \c i2. Index -1 indicates the last element.
    const Vec<Num_T> operator()(int i1, int i2) const;
    //! Sub-vector where the elements are given by the list \c indexlist
    const Vec<Num_T> operator()(const Vec<int> &indexlist) const;

    //! Accessor-style method. First element is 0
    const Num_T &get(int i) const { it_assert0(i>=0&&i<datasize, "method get()"); return data[i]; }
    //! Sub-vector with elements from \c i1 to \c i2. Index -1 indicates the last element.
    const Vec<Num_T> get(int i1, int i2) const;
    //! Modifier-style method. First element is 0
    void set(int i, const Num_T &v) { it_assert0(i>=0&&i<datasize, "method set()"); data[i]=v; }

    //! Matrix transpose. Converts to a matrix with the vector in the first row
    Mat<Num_T> transpose() const;
    //! Matrix transpose. Converts to a matrix with the vector in the first row
    Mat<Num_T> T() const { return this->transpose(); }
    //! Hermitian matrix transpose. Converts to a matrix with the conjugate of the vector in the first row
    Mat<Num_T> hermitian_transpose() const;
    //! Hermitian matrix transpose. Converts to a matrix with the conjugate of the vector in the first row
    Mat<Num_T> H() const { return this->hermitian_transpose(); }

    //! Addition of vector
    void operator+=(const Vec<Num_T> &v);
    //! Addition of scalar
    void operator+=(Num_T t) { for (int i=0;i<datasize;i++) data[i]+=t; }
    //! Addition of two vectors
    friend Vec<Num_T> operator+<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Addition of a vector and a scalar
    friend Vec<Num_T> operator+<>(const Vec<Num_T> &v, Num_T t);
    //! Addition of a scalar and a vector
    friend Vec<Num_T> operator+<>(Num_T t, const Vec<Num_T> &v);

    //! Subtraction of vector
    void operator-=(const Vec<Num_T> &v);
    //! Subtraction of scalar
    void operator-=(Num_T t) { for (int i=0;i<datasize;i++) data[i]-=t; }
    //! Subtraction of \c v2 from \c v1
    friend Vec<Num_T> operator-<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Subtraction of scalar from vector
    friend Vec<Num_T> operator-<>(const Vec<Num_T> &v, Num_T t);
    //! Sutraction of vector from scalar
    friend Vec<Num_T> operator-<>(Num_T t, const Vec<Num_T> &v);
    //! Negation of vector
    friend Vec<Num_T> operator-<>(const Vec<Num_T> &v);

    //! Multiply with a scalar
    void operator*=(Num_T t) { for (int i=0;i<datasize;i++) data[i] *= t; }
    //! Inner (dot) product
    friend Num_T operator*<>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Inner (dot) product
    friend Num_T dot <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Outer product of two vectors v1 and v2
    friend Mat<Num_T> outer_product <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Elementwise multiplication of vector and scalar
    friend Vec<Num_T> operator*<>(const Vec<Num_T> &v, Num_T t);
    //! Elementwise multiplication of vector and scalar
    friend Vec<Num_T> operator*<>(Num_T t, const Vec<Num_T> &v);
    //! Elementwise multiplication
    friend Vec<Num_T> elem_mult <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Elementwise multiplication of three vectors
    friend Vec<Num_T> elem_mult <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3);
    //! Elementwise multiplication of four vectors
    friend Vec<Num_T> elem_mult <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3, const Vec<Num_T> &v4);

    //! Elementwise division
    void operator/=(Num_T t) { for (int i=0;i<datasize;i++) data[i]/=t; }
    //! Elementwise division
    friend Vec<Num_T> operator/<>(const Vec<Num_T> &v, Num_T t);
    //! Elementwise division
    friend Vec<Num_T> operator/<>(const Num_T t, const Vec<Num_T> &v);
    //! Elementwise division
    void operator/=(const Vec<Num_T> &v);
    //! Elementwise division
    friend Vec<Num_T> elem_div <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2);
    //! Elementwise division
    friend Vec<Num_T> elem_div <>(const Num_T t, const Vec<Num_T> &v);

    //! Get the elements in the vector where \c binlist is \c 1
    Vec<Num_T> get(const Vec<bin> &binlist) const;
    //! Get the right \c nr elements from the vector
    Vec<Num_T> right(int nr) const;
    //! Get the left \c nr elements from the vector
    Vec<Num_T> left(int nr) const;
    //! Get the middle part of vector from \c start including \c nr elements
    Vec<Num_T> mid(int start, int nr) const;
    //! Split the vector into two parts at element \c pos. Return the first part and keep the second.
    Vec<Num_T> split(int pos);
    //! Shift in element \c In at position 0 \c n times
    void shift_right(Num_T In, int n=1);
    //! Shift in vector \c In at position 0
    void shift_right(const Vec<Num_T> &In);
    //! Shift out the \c n left elements and a the same time shift in the element \c at last position \c n times
    void shift_left(Num_T In, int n=1);
    //! Shift in vector \c In at last position
    void shift_left(const Vec<Num_T> &In);

    //! Append element \c a to the end of the vector \c v
    friend Vec<Num_T> concat <>(const Vec<Num_T> &v, const Num_T a);
    //! Concat element \c a to the beginning of the vector \c v
    friend Vec<Num_T> concat <>(const Num_T a, const Vec<Num_T> &v);
    //! Concat vectors \c v1 and \c v2
    friend Vec<Num_T> concat <>(const Vec<Num_T> &v1,const Vec<Num_T> &v2);
    //! Concat vectors \c v1, \c v2 and \c v3
    friend Vec<Num_T> concat <>(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3);

    //! Set subvector defined by indicies \c i1 to \c i2 to vector \c v
    void set_subvector(int i1, int i2, const Vec<Num_T> &v);
    //! Set subvector defined by first index \c i and size of vector \c v to \c v
    void set_subvector(int i, const Vec<Num_T> &v);
    //! Set subvector defined by indicies i1 to i2 to constant t
    void set_subvector(int i1, int i2, const Num_T t);
    //! Replace the elements from \c pos by the vector \c v
    void replace_mid(int pos, const Vec<Num_T> &v);
    //! Delete element number \c index
    void del(int index);
    //! Delete elements from \c i1 to \c i2
    void del(int i1, int i2);
    //! Insert element \c in at \c index
    void ins(int index, Num_T in);
    //! Insert vector \c in at \c index
    void ins(int index, const Vec<Num_T> &in);

    //! Assign all elements in vector to \c t
    void operator=(Num_T t) { for (int i=0;i<datasize;i++) data[i] = t; }
    //! Assign vector the value and length of \c v
    void operator=(const Vec<Num_T> &v);
    //! Assign vector equal to the 1-dimensional matrix \c m
    void operator=(const Mat<Num_T> &m);
    //! Assign vector the values in the string \c values
    void operator=(const char *values) { set(values); }

    //! Elementwise equal to the scalar
    Vec<bin> operator==(const Num_T value) const;
    //! Elementwise not-equal to the scalar
    Vec<bin> operator!=(const Num_T value) const;
    //! Elementwise less than the scalar
    Vec<bin> operator<(const Num_T value) const;
    //! Elementwise less than and equal to the scalar
    Vec<bin> operator<=(const Num_T value) const;
    //! Elementwise greater than the scalar
    Vec<bin> operator>(const Num_T value) const;
    //! Elementwise greater than and equal to the scalar
    Vec<bin> operator>=(const Num_T value) const;

    //! Compare two vectors. False if wrong sizes or different values
    bool operator==(const Vec<Num_T> &v) const;
    //! Compare two vectors. True if different
    bool operator!=(const Vec<Num_T> &v) const;

    //! Index operator without boundary check. Not recommended to use.
    Num_T &_elem(int i) { return data[i]; }
    //! Index operator without boundary check. Not recommended to use.
    const Num_T &_elem(int i) const { return data[i]; }

    //! Get the pointer to the internal structure. Not recommended to use.
    Num_T *_data() { return data; }

    //! Get the pointer to the internal structure. Not recommended to use.
    const Num_T *_data() const { return data; }

  protected:
    //! Allocate storage for a vector of length \c size.
    void alloc(int size)
    {
      if ( datasize == size ) return;

      free();  // Free memory (if any allocated)
      if (size == 0) return;

      create_elements(data, size, factory);
      datasize=size;
      it_assert1(data, "Out of memory in Vec::alloc()");
    }

    //! Free the storage space allocated by the vector
    void free() { delete [] data;  init(); }

    //! The current number of elements in the vector
    int datasize;
    //! A pointer to the data area
    Num_T *data;
    //! Element factory (set to DEFAULT_FACTORY to use Num_T default constructors only)
    const Factory &factory;

  private:
    void init() { data = 0; datasize = 0; }
  };

  //-----------------------------------------------------------------------------------
  // Type definitions of vec, cvec, ivec, svec, and bvec
  //-----------------------------------------------------------------------------------

  /*!
    \relates Vec
    \brief Definition of double vector type
  */
  typedef Vec<double> vec;

  /*!
    \relates Vec
    \brief Definition of complex<double> vector type
  */
  typedef Vec<std::complex<double> > cvec;

  /*!
    \relates Vec
    \brief Definition of integer vector type
  */
  typedef Vec<int> ivec;

  /*!
    \relates Vec
    \brief Definition of short vector type
  */
  typedef Vec<short int> svec;

  /*!
    \relates Vec
    \brief Definition of binary vector type
  */
  typedef Vec<bin> bvec;

} //namespace itpp

#include "itpp/base/mat.h"

namespace itpp {

  //-----------------------------------------------------------------------------------
  // Declaration of input and output streams for Vec
  //-----------------------------------------------------------------------------------

  /*!
    \relates Vec
    \brief Stream output of vector
  */
  template <class Num_T>
    std::ostream &operator<<(std::ostream &os, const Vec<Num_T> &v);

  /*!
    \relates Vec
    \brief Stream input of vector

    The input can be on the form "1 2 3" or "[1 2 3]", i.e. with or without
    brackets. The first form is compatible with the set method, while the
    second form is compatible with the ostream operator. The elements can be
    separated by blank space or commas. "[]" means an empty vector. "1:4"
    means "1 2 3 4". "1:3:10" means every third integer from 1 to 10, i.e.
    "1 4 7 10".
  */
  template <class Num_T>
    std::istream &operator>>(std::istream &is, Vec<Num_T> &v);

  //-----------------------------------------------------------------------------------
  // Implementation of templated Vec members and friends
  //-----------------------------------------------------------------------------------

  template<class Num_T> inline
    Vec<Num_T>::Vec(const Vec<Num_T> &v) : factory(v.factory)
  {
    init();
    alloc(v.datasize);
    copy_vector(datasize, v.data, data);
  }

  template<class Num_T> inline
    Vec<Num_T>::Vec(const Vec<Num_T> &v, const Factory &f) : factory(f)
  {
    init();
    alloc(v.datasize);
    copy_vector(datasize, v.data, data);
  }

  template<class Num_T>
    void Vec<Num_T>::set_size(int size, bool copy)
  {
    it_assert1(size >= 0, "New size must not be negative in Vec::set_size()");
    if (size!=datasize) {
      if (copy) {
	Vec<Num_T> temp(*this);

	alloc(size);
	for (int i=0; i<size; i++)
	  data[i] = i < temp.datasize ? temp.data[i] : Num_T(0);
      } else
	alloc(size);
    }
  }

  template<> bool cvec::set(const char *values);
  template<> bool bvec::set(const char *values);

  template<class Num_T>
    bool Vec<Num_T>::set(const char *values)
  {
    std::istringstream buffer(values);
    Num_T b, c;
    int pos=0, maxpos=10;
    
    alloc(maxpos);
    
    while (buffer.peek()!=EOF) {

      switch (buffer.peek()) {
      case ':': // reads format a:b:c or a:b
	buffer.get();
	if (!buffer.eof()) {
	  buffer >> b;
	}
	if (!buffer.eof() && buffer.peek() == ':') {
	  buffer.get();
	  if (!buffer.eof()) {
	    buffer >> c;

	    while (sign(b)*(data[pos-1]+b-c)<=0) {
	      pos++;
	      if (pos > maxpos) {
		maxpos=maxpos*2;
		set_size(maxpos, true);
	      }
	      data[pos-1]=data[pos-2]+b;
	    }
	  }
	} else {
	  while (data[pos-1]<b) {
	    pos++;
	    if (pos > maxpos) {
	      maxpos=maxpos*2;
	      set_size(maxpos, true);
	    }
	    data[pos-1]=data[pos-2]+1;
	  }
	}
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

  template<class Num_T>
    bool Vec<Num_T>::set(const std::string &str)
  {
    return set(str.c_str());
  }

  template<class Num_T> inline
    const Vec<Num_T> Vec<Num_T>::operator()(int i1, int i2) const
  {
    if (i1 == -1)	i1 = datasize-1;
    if (i2 == -1) i2 = datasize-1;

    it_assert1(i1>=0 && i2>=0 && i1<datasize && i2<datasize, "Vec<Num_T>::operator()(i1,i2): indicies out of range");
    it_assert1(i2>=i1, "Vec<Num_T>::op(i1,i2): i2 >= i1 necessary");

    Vec<Num_T> s(i2-i1+1);
    copy_vector(s.datasize, data+i1, s.data);

    return s;
  }

  template<class Num_T> inline
    const Vec<Num_T> Vec<Num_T>::get(int i1, int i2) const
  {
    return (*this)(i1, i2);
  }

  template<class Num_T>
    const Vec<Num_T> Vec<Num_T>::operator()(const Vec<int> &indexlist) const
  {
    Vec<Num_T> temp(indexlist.length());
    for (int i=0;i<indexlist.length();i++) {
      it_assert((indexlist(i)>=0) && (indexlist(i) < datasize), "Vec<Num_T>::operator()(ivec &): index outside range");
      temp(i)=data[indexlist(i)];
    }
    return temp;
  }

  template<class Num_T>
    Mat<Num_T> Vec<Num_T>::transpose() const
  {
    Mat<Num_T> temp(1, datasize);
    for (int i=0; i<datasize; i++)
      temp(i) = data[i];

    return temp;
  }

  template<> Mat<std::complex<double> > cvec::hermitian_transpose() const;

  template<class Num_T>
    Mat<Num_T> Vec<Num_T>::hermitian_transpose() const
  {
    Mat<Num_T> temp(1, datasize);
    for (int i=0; i<datasize; i++)
      temp(i) = data[i];

    return temp;
  }

  template<class Num_T> inline
    void Vec<Num_T>::operator+=(const Vec<Num_T> &v)
  {
    int i;

    if (datasize == 0) { // if not assigned a size.
      alloc(v.datasize);
      for (i=0; i<v.datasize; i++)
	data[i] = v.data[i];
    } else {
      it_assert1(datasize==v.datasize, "Vec<Num_T>::operator+=: wrong sizes");
      for (i=0; i<datasize; i++)
	data[i] += v.data[i];
    }
  }

  template<class Num_T> inline
    Vec<Num_T> operator+(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Vec<Num_T> r(v1.datasize);

    it_assert1(v1.datasize==v2.datasize, "Vec<Num_T>::operator+: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r.data[i] = v1.data[i] + v2.data[i];

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator+(const Vec<Num_T> &v, Num_T t)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = v.data[i] + t;

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator+(Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t + v.data[i];

    return r;
  }

  template<class Num_T> inline
    void Vec<Num_T>::operator-=(const Vec<Num_T> &v)
  {
    int i;

    if (datasize == 0) { // if not assigned a size.
      alloc(v.datasize);
      for (i=0; i<v.datasize; i++)
	data[i] = -v.data[i];
    } else {
      it_assert1(datasize==v.datasize, "Vec<Num_T>::operator-=: wrong sizes");
      for (i=0; i<datasize; i++)
	data[i] -= v.data[i];
    }
  }

  template<class Num_T> inline
    Vec<Num_T> operator-(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Vec<Num_T> r(v1.datasize);

    it_assert1(v1.datasize==v2.datasize, "Vec<Num_T>::operator-: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r.data[i] = v1.data[i] - v2.data[i];

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator-(const Vec<Num_T> &v, Num_T t)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = v.data[i] - t;

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator-(Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t - v.data[i];

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator-(const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = -v.data[i];

    return r;
  }

#if defined(HAVE_CBLAS) || defined(HAVE_MKL)
 template<> double dot(const Vec<double> &v1, const Vec<double> &v2);
 template<> std::complex<double> dot(const Vec<std::complex<double> > &v1, const Vec<std::complex<double> > &v2);
#endif

  template<class Num_T> inline
    Num_T dot(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Num_T r=Num_T(0);

    it_assert1(v1.datasize==v2.datasize, "Vec<Num_T>::dot: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r += v1.data[i] * v2.data[i];

    return r;
  }

  template<class Num_T> inline
    Mat<Num_T> outer_product(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i, j;

    it_assert1(v1.datasize>0 && v2.datasize>0, "outer_product:: Vector of zero size");

    Mat<Num_T> r(v1.datasize, v2.datasize);

    for (i=0; i<v1.datasize; i++) {
      for (j=0; j<v2.datasize; j++) {
	r(i,j) = v1.data[i] * v2.data[j];
      }
    }

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator*(const Vec<Num_T> &v, Num_T t)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = v.data[i] * t;

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator*(Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t * v.data[i];

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> elem_mult(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Vec<Num_T> r(v1.datasize);

    it_assert1(v1.datasize==v2.datasize, "Vec<Num_T>::elem_mult: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r.data[i] = v1.data[i] * v2.data[i];

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> elem_mult(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3)
  {
    int i;
    Vec<Num_T> r(v1.datasize);

    it_assert1(v1.datasize==v2.datasize, "Vec<Num_T>::elem_mult: wrong sizes");
    it_assert1(v2.datasize==v3.datasize, "Vec<Num_T>::elem_mult: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r.data[i] = v1.data[i] * v2.data[i] * v3.data[i];

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> elem_mult(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3, const Vec<Num_T> &v4)
  {
    int i;
    Vec<Num_T> r(v1.datasize);

    it_assert1(v1.datasize==v2.datasize, "Vec<Num_T>::elem_mult: wrong sizes");
    it_assert1(v2.datasize==v3.datasize, "Vec<Num_T>::elem_mult: wrong sizes");
    it_assert1(v3.datasize==v4.datasize, "Vec<Num_T>::elem_mult: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r.data[i] = v1.data[i] * v2.data[i] * v3.data[i] * v4.data[i];

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator/(const Vec<Num_T> &v, Num_T t)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = v.data[i] / t;

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> operator/(const Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t / v.data[i];

    return r;
  }

  template<class Num_T> inline
    void Vec<Num_T>::operator/=(const Vec<Num_T> &v)
  {
    int i;

    it_assert1(datasize==v.datasize, "Vec<Num_T>::operator/=: wrong sizes");
    for (i=0; i<datasize; i++)
      data[i] /= v.data[i];
  }

  template<class Num_T> inline
    Vec<Num_T> elem_div(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Vec<Num_T> r(v1.datasize);

    it_assert1(v1.datasize==v2.datasize, "elem_div: wrong sizes");
    for (i=0; i<v1.datasize; i++)
      r.data[i] = v1.data[i] / v2.data[i];

    return r;
  }

  template<class Num_T> inline
    Vec<Num_T> elem_div(const Num_T t, const Vec<Num_T> &v)
  {
    int i;
    Vec<Num_T> r(v.datasize);

    for (i=0; i<v.datasize; i++)
      r.data[i] = t / v.data[i];

    return r;
  }

  template<class Num_T>
    Vec<Num_T> Vec<Num_T>::get(const Vec<bin> &binlist) const
  {
    it_assert1(datasize == binlist.size(), "Vec<Num_T>::get(bvec &): wrong sizes");
    Vec<Num_T> temp(binlist.length());
    int j=0;

    for (int i=0;i<binlist.length();i++) {
      if (binlist(i) == bin(1)) {
	temp(j)=data[i];
	j++;
      }
    }
    temp.set_size(j, true);
    return temp;
  }

  template<class Num_T> inline
    Vec<Num_T> Vec<Num_T>::right(int nr) const
  {
    it_assert1(nr<=datasize, "Vec<Num_T>::right: index out of range");
    Vec<Num_T> temp(nr);
    if (nr!=0) {
      copy_vector(nr, &data[datasize-nr], &temp[0]);
    }
    return temp;
  }

  template<class Num_T> inline
    Vec<Num_T> Vec<Num_T>::left(int nr) const
  {
    it_assert1(nr<=datasize, "Vec<Num_T>::left: index out of range");
    Vec<Num_T> temp(nr);
    if (nr!=0) {
      copy_vector(nr, &data[0], &temp[0]);
    }
    return temp;
  }

  template<class Num_T> inline
    Vec<Num_T> Vec<Num_T>::mid(int start, int nr) const
  {
    it_assert1((start>=0)&& ((start+nr)<=datasize), "Vec<Num_T>::mid: indexing out of range");
    Vec<Num_T> temp(nr);

    if (nr!=0) {
      copy_vector(nr, &data[start], &temp[0]);
    }
    return temp;
  }

  template<class Num_T>
    Vec<Num_T> Vec<Num_T>::split(int Position)
  {
    it_assert1((Position>=0) && (Position<=datasize), "Vec<Num_T>::split: index out of range");
    Vec<Num_T> Temp1(Position);
    Vec<Num_T> Temp2(datasize-Position);
    int	 i;

    for (i=0;i<Position;i++) {
      Temp1[i]=data[i];
    }
    for (i=Position;i<datasize;i++) {
      Temp2[i-Position]=data[i];
    }
    (*this)=Temp2;
    return Temp1;
  }

  template<class Num_T>
    void Vec<Num_T>::shift_right(Num_T In, int n)
  {
    int i=datasize;

    it_assert1(n>=0, "Vec<Num_T>::shift_right: index out of range");
    while (--i >= n)
      data[i] = data[i-n];
    while (i >= 0)
      data[i--] = In;
  }

  template<class Num_T>
    void Vec<Num_T>::shift_right(const Vec<Num_T> &In)
  {
    int	i;

    for (i=datasize-1; i>=In.datasize; i--)
      data[i]=data[i-In.datasize];
    for (i=0; i<In.datasize; i++)
      data[i]=In[i];
  }

  template<class Num_T>
    void Vec<Num_T>::shift_left(Num_T In, int n)
  {
    int i;

    it_assert1(n>=0, "Vec<Num_T>::shift_left: index out of range");
    for (i=0; i<datasize-n; i++)
      data[i] = data[i+n];
    while (i < datasize)
      data[i++] = In;
  }

  template<class Num_T>
    void Vec<Num_T>::shift_left(const Vec<Num_T> &In)
  {
    int	i;

    for (i=0; i<datasize-In.datasize; i++)
      data[i]=data[i+In.datasize];
    for (i=datasize-In.datasize; i<datasize; i++)
      data[i]=In[i-datasize+In.datasize];
  }

  template<class Num_T>
    Vec<Num_T> concat(const Vec<Num_T> &v, const Num_T a)
  {
    Vec<Num_T> temp(v.size()+1);

    for (int i=0; i<v.size(); i++)
      temp(i) = v(i);
    temp(v.size()) = a;

    return temp;
  }

  template<class Num_T>
    Vec<Num_T> concat(const Num_T a, const Vec<Num_T> &v)
  {
    Vec<Num_T> temp(v.size()+1);

    temp(0) = a;

    for (int i=0; i<v.size(); i++)
      temp(i+1) = v(i);

    return temp;
  }

  template<class Num_T>
    Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2)
  {
    int i;
    Vec<Num_T> temp(v1.size()+v2.size());

    for (i=0;i<v1.size();i++) {
      temp[i] = v1[i];
    }
    for (i=0;i<v2.size();i++) {
      temp[v1.size()+i] = v2[i];
    }
    return temp;
  }

  template<class Num_T>
    Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3)
  {
    // There should be some error control?
    int i;
    Vec<Num_T> temp(v1.size()+v2.size()+v3.size());

    for (i=0;i<v1.size();i++) {
      temp[i] = v1[i];
    }
    for (i=0;i<v2.size();i++) {
      temp[v1.size()+i] = v2[i];
    }
    for (i=0;i<v3.size();i++) {
      temp[v1.size()+v2.size()+i] = v3[i];
    }
    return temp;
  }

  template<class Num_T>
    Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3, const Vec<Num_T> &v4)
  {
    // There should be some error control?
    int i;
    Vec<Num_T> temp(v1.size()+v2.size()+v3.size()+v4.size());

    for (i=0;i<v1.size();i++) {
      temp[i] = v1[i];
    }
    for (i=0;i<v2.size();i++) {
      temp[v1.size()+i] = v2[i];
    }
    for (i=0;i<v3.size();i++) {
      temp[v1.size()+v2.size()+i] = v3[i];
    }
    for (i=0;i<v4.size();i++) {
      temp[v1.size()+v2.size()+v3.size()+i] = v4[i];
    }
    return temp;
  }

  template<class Num_T>
    Vec<Num_T> concat(const Vec<Num_T> &v1, const Vec<Num_T> &v2, const Vec<Num_T> &v3, const Vec<Num_T> &v4, const Vec<Num_T> &v5)
  {
    // There should be some error control?
    int i;
    Vec<Num_T> temp(v1.size()+v2.size()+v3.size()+v4.size()+v5.size());

    for (i=0;i<v1.size();i++) {
      temp[i] = v1[i];
    }
    for (i=0;i<v2.size();i++) {
      temp[v1.size()+i] = v2[i];
    }
    for (i=0;i<v3.size();i++) {
      temp[v1.size()+v2.size()+i] = v3[i];
    }
    for (i=0;i<v4.size();i++) {
      temp[v1.size()+v2.size()+v3.size()+i] = v4[i];
    }
    for (i=0;i<v5.size();i++) {
      temp[v1.size()+v2.size()+v3.size()+v4.size()+i] = v5[i];
    }
    return temp;
  }

  template<class Num_T> inline
    void Vec<Num_T>::set_subvector(int i1, int i2, const Vec<Num_T> &v)
  {
    if (i1 == -1) i1 = datasize-1;
    if (i2 == -1) i2 = datasize-1;

    it_assert1(i1>=0 && i2>=0 && i1<datasize && i2<datasize, "Vec<Num_T>::set_subvector(): indicies out of range");
    it_assert1(i2>=i1, "Vec<Num_T>::set_subvector(): i2 >= i1 necessary");
    it_assert1(i2-i1+1 == v.datasize, "Vec<Num_T>::set_subvector(): wrong sizes");

    copy_vector(v.datasize, v.data, data+i1);
  }

  template<class Num_T> inline
    void Vec<Num_T>:: set_subvector(int i, const Vec<Num_T> &v)
  {
    it_assert1(i>=0, "Vec<Num_T>::set_subvector(): index out of range");
    it_assert1(i+v.datasize <= datasize, "Vec<Num_T>::set_subvector(): too long input vector");
    copy_vector(v.datasize, v.data, data+i);
  }

  template<class Num_T>
    void Vec<Num_T>::set_subvector(int i1, int i2, const Num_T t)
  {
    if (i1 == -1) i1 = datasize-1;
    if (i2 == -1) i2 = datasize-1;

    it_assert1(i1>=0 && i2>=0 && i1<datasize && i2<datasize, "Vec<Num_T>::set_subvector(): indicies out of range");
    it_assert1(i2>=i1, "Vec<Num_T>::set_subvector(): i2 >= i1 necessary");

    for (int i=i1;i<=i2;i++)
      data[i] = t;
  }

  template<class Num_T>
    void Vec<Num_T>::replace_mid(int pos, const Vec<Num_T> &v)
  {
    it_assert1((pos>=0) && ((pos+v.length())<=datasize), "Vec<Num_T>::replace_mid: indexing out of range");
    copy_vector(v.datasize, v.data, &data[pos]);
  }

  template<class Num_T>
    void Vec<Num_T>::del(int index)
  {
    it_assert1((index>=0) && (index<datasize), "Vec<Num_T>::del: index out of range");
    Vec<Num_T> Temp(*this);
    int i;

    set_size(datasize-1, false);
    for (i=0;i<index;i++) {
      data[i]=Temp[i];
    }
    for (i=index;i<datasize;i++) {
      data[i]=Temp[i+1];
    }
  }

  template<class Num_T>
    void  Vec<Num_T>::del(int i1, int i2) 
  {
    it_assert1((i1>=0) && (i2<datasize) && (i1<i2), "Vec<Num_T>::del: index out of range");

    Vec<Num_T> Temp(*this);
    int new_size = datasize-(i2-i1+1);
    set_size(new_size, false);
    copy_vector(i1, Temp.data, data);
    copy_vector(datasize-i1, &Temp.data[i2+1], &data[i1]);
  }

  template<class Num_T>
    void Vec<Num_T>::ins(int index, Num_T in)
  {
    it_assert1((index>=0) && (index<=datasize), "Vec<Num_T>::ins: index out of range");
    Vec<Num_T> Temp(*this);

    set_size(datasize+1, false);
    copy_vector(index, Temp.data, data);
    data[index]=in;
    copy_vector(Temp.datasize-index, Temp.data+index, data+index+1);
  }

  template<class Num_T>
    void Vec<Num_T>::ins(int index, const Vec<Num_T> &in)
  {
    it_assert1((index>=0) && (index<=datasize), "Vec<Num_T>::ins: index out of range");
    Vec<Num_T> Temp(*this);

    set_size(datasize+in.length(), false);
    copy_vector(index, Temp.data, data);
    copy_vector(in.size(), in.data, &data[index]);
    copy_vector(Temp.datasize-index, Temp.data+index, data+index+in.size());
  }

  template<class Num_T> inline
    void Vec<Num_T>::operator=(const Vec<Num_T> &v)
    {
      set_size(v.datasize, false);
      copy_vector(datasize, v.data, data);
    }

  template<class Num_T> inline
    void Vec<Num_T>::operator=(const Mat<Num_T> &m)
    {
      it_assert1( (m.cols() == 1 && datasize == m.rows()) ||
		  (m.rows() == 1 && datasize == m.cols()), "vec::op=(mat); wrong size");

      if (m.cols() == 1) {
	set_size(m.rows(), false);
	copy_vector(m.rows(), m._data(), data);
      } else if (m.rows() == 1) {
	set_size(m.cols(), false);
	copy_vector(m.cols(), m._data(), m.rows(), data, 1);
      } else
	it_error("vec::op=(mat); wrong size");
    }

  template<> bvec cvec::operator==(const std::complex<double>) const;

  template<class Num_T>
    bvec Vec<Num_T>::operator==(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec<Num_T>::operator==: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)==value);

    return temp;
  }

  template<> bvec cvec::operator!=(const std::complex<double>) const;

  template<class Num_T>
    bvec Vec<Num_T>::operator!=(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec<Num_T>::operator!=: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)!=value);

    return temp;
  }

  template<> bvec cvec::operator<(const std::complex<double>) const;

  template<class Num_T>
    bvec Vec<Num_T>::operator<(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec<Num_T>::operator<: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)<value);

    return temp;
  }

  template<> bvec cvec::operator<=(const std::complex<double>) const;

  template<class Num_T>
    bvec Vec<Num_T>::operator<=(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec<Num_T>::operator<=: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)<=value);

    return temp;
  }

  template<> bvec cvec::operator>(const std::complex<double>) const;

  template<class Num_T>
    bvec Vec<Num_T>::operator>(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec<Num_T>::operator>: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)>value);

    return temp;
  }

  template<> bvec cvec::operator>=(const std::complex<double>) const;

  template<class Num_T>
    bvec Vec<Num_T>::operator>=(const Num_T value) const
  {
    it_assert(datasize > 0, "Vec<Num_T>::operator>=: vector must have size > 0");
    Vec<Num_T> invector(*this);
    bvec temp(invector.length());

    for (int i=0;i<invector.length();i++)
      temp(i)=(invector(i)>=value);

    return temp;
  }

  template<class Num_T>
    bool Vec<Num_T>::operator==(const Vec<Num_T> &invector) const
  {
    // OBS ! if wrong size, return false
    if (datasize!=invector.datasize) return false;
    for (int i=0;i<datasize;i++) {
      if (data[i]!=invector.data[i]) return false;
    }
    return true;
  }

  template<class Num_T>
    bool Vec<Num_T>::operator!=(const Vec<Num_T> &invector) const
  {
    if (datasize!=invector.datasize) return true;
    for (int i=0;i<datasize;i++) {
      if (data[i]!=invector.data[i]) return true;
    }
    return false;
  }

  template <class Num_T>
    std::ostream &operator<<(std::ostream &os, const Vec<Num_T> &v)
  {
    int i, sz=v.length();

    os << "[" ;
    for (i=0; i<sz; i++) {
      os << v(i) ;
      if (i < sz-1)
	os << " ";
    }
    os << "]" ;

    return os;
  }

  template <class Num_T>
    std::istream &operator>>(std::istream &is, Vec<Num_T> &v)
  {
    std::ostringstream buffer;
    bool started = false;
    bool finished = false;
    bool brackets = false;
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
          if (started) {
            // Unexpected left bracket
            is.setstate(std::ios_base::failbit);
            finished = true;
          } else {
            started = true;
            brackets = true;
          }
        } else if (c == ']') {
          if (!started || !brackets) {
            // Unexpected right bracket
            is.setstate(std::ios_base::failbit);
            finished = true;
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
      v.set_size(0, false);
    } else {
      v.set(buffer.str());
    }

    return is;
  }

  //---------------------------------------------------------------------
  // Instantiations
  //---------------------------------------------------------------------

  //--------- class instantiations -------------

#ifndef _MSC_VER
  //! Template instantiation of Vec<double>
  extern template class Vec<double>;
  //! Template instantiation of Vec<int>
  extern template class Vec<int>;
  //! Template instantiation of Vec<short int>
  extern template class Vec<short int>;
  //! Template instantiation of Vec<complex<double> >
  extern template class Vec<std::complex<double> >;
  //! Template instantiation of Vec<bin>
  extern template class Vec<bin>;

  //------------- Addition operator ----------

  //! Template instantiation of operator+
  extern template vec operator+(const vec &v1, const vec &v2);
  //! Template instantiation of operator+
  extern template cvec operator+(const cvec &v1, const cvec &v2);
  //! Template instantiation of operator+
  extern template ivec operator+(const ivec &v1, const ivec &v2);
  //! Template instantiation of operator+
  extern template svec operator+(const svec &v1, const svec &v2);
  //! Template instantiation of operator+
  extern template bvec operator+(const bvec &v1, const bvec &v2);

  //! Template instantiation of operator+
  extern template vec operator+(const vec &v1, double t);
  //! Template instantiation of operator+
  extern template cvec operator+(const cvec &v1, std::complex<double> t);
  //! Template instantiation of operator+
  extern template ivec operator+(const ivec &v1, int t);
  //! Template instantiation of operator+
  extern template svec operator+(const svec &v1, short t);
  //! Template instantiation of operator+
  extern template bvec operator+(const bvec &v1, bin t);

  //! Template instantiation of operator+
  extern template vec operator+(double t, const vec &v1);
  //! Template instantiation of operator+
  extern template cvec operator+(std::complex<double> t, const cvec &v1);
  //! Template instantiation of operator+
  extern template ivec operator+(int t, const ivec &v1);
  //! Template instantiation of operator+
  extern template svec operator+(short t, const svec &v1);
  //! Template instantiation of operator+
  extern template bvec operator+(bin t, const bvec &v1);

  //------------- Subraction operator ----------

  //! Template instantiation of operator-
  extern template vec operator-(const vec &v1, const vec &v2);
  //! Template instantiation of operator-
  extern template cvec operator-(const cvec &v1, const cvec &v2);
  //! Template instantiation of operator-
  extern template ivec operator-(const ivec &v1, const ivec &v2);
  //! Template instantiation of operator-
  extern template svec operator-(const svec &v1, const svec &v2);
  //! Template instantiation of operator-
  extern template bvec operator-(const bvec &v1, const bvec &v2);

  //! Template instantiation of operator-
  extern template vec operator-(const vec &v, double t);
  //! Template instantiation of operator-
  extern template cvec operator-(const cvec &v, std::complex<double> t);
  //! Template instantiation of operator-
  extern template ivec operator-(const ivec &v, int t);
  //! Template instantiation of operator-
  extern template svec operator-(const svec &v, short t);
  //! Template instantiation of operator-
  extern template bvec operator-(const bvec &v, bin t);

  //! Template instantiation of operator-
  extern template vec operator-(double t, const vec &v);
  //! Template instantiation of operator-
  extern template cvec operator-(std::complex<double> t, const cvec &v);
  //! Template instantiation of operator-
  extern template ivec operator-(int t, const ivec &v);
  //! Template instantiation of operator-
  extern template svec operator-(short t, const svec &v);
  //! Template instantiation of operator-
  extern template bvec operator-(bin t, const bvec &v);

  //---------- Unary minus -------------

  //! Template instantiation of operator-
  extern template vec operator-(const vec &v);
  //! Template instantiation of operator-
  extern template cvec operator-(const cvec &v);
  //! Template instantiation of operator-
  extern template ivec operator-(const ivec &v);
  //! Template instantiation of operator-
  extern template svec operator-(const svec &v);
  //! Template instantiation of operator-
  extern template bvec operator-(const bvec &v);

  //------------- Multiplication operator ----------

  //! Template instantiation of dot
  extern template double dot(const vec &v1, const vec &v2);
  //! Template instantiation of dot
  extern template std::complex<double> dot(const cvec &v1, const cvec &v2);
  //! Template instantiation of dot
  extern template int dot(const ivec &v1, const ivec &v2);
  //! Template instantiation of dot
  extern template short dot(const svec &v1, const svec &v2);
  //! Template instantiation of dot
  extern template bin dot(const bvec &v1, const bvec &v2);

  //! Template instantiation of operator*
  extern template int operator*(const ivec &v1, const ivec &v2);
  //! Template instantiation of operator*
  extern template short operator*(const svec &v1, const svec &v2);
  //! Template instantiation of operator*
  extern template bin operator*(const bvec &v1, const bvec &v2);

  //! Template instantiation of outer_product
  extern template mat outer_product(const vec &v1, const vec &v2);
  //! Template instantiation of outer_product
  extern template cmat outer_product(const cvec &v1, const cvec &v2);
  //! Template instantiation of outer_product
  extern template imat outer_product(const ivec &v1, const ivec &v2);
  //! Template instantiation of outer_product
  extern template smat outer_product(const svec &v1, const svec &v2);
  //! Template instantiation of outer_product
  extern template bmat outer_product(const bvec &v1, const bvec &v2);

  //! Template instantiation of operator*
  extern template vec operator*(const vec &v, double t);
  //! Template instantiation of operator*
  extern template cvec operator*(const cvec &v, std::complex<double> t);
  //! Template instantiation of operator*
  extern template ivec operator*(const ivec &v, int t);
  //! Template instantiation of operator*
  extern template svec operator*(const svec &v, short t);
  //! Template instantiation of operator*
  extern template bvec operator*(const bvec &v, bin t);

  //! Template instantiation of operator*
  extern template vec operator*(double t, const vec &v);
  //! Template instantiation of operator*
  extern template cvec operator*(std::complex<double> t, const cvec &v);
  //! Template instantiation of operator*
  extern template ivec operator*(int t, const ivec &v);
  //! Template instantiation of operator*
  extern template svec operator*(short t, const svec &v);
  //! Template instantiation of operator*
  extern template bvec operator*(bin t, const bvec &v);

  //------------- Elementwise Multiplication operator (two vectors) ----------

  //! Template instantiation of elem_mult
  extern template vec elem_mult(const vec &v1, const vec &v2);
  //! Template instantiation of elem_mult
  extern template cvec elem_mult(const cvec &v1, const cvec &v2);
  //! Template instantiation of elem_mult
  extern template ivec elem_mult(const ivec &v1, const ivec &v2);
  //! Template instantiation of elem_mult
  extern template svec elem_mult(const svec &v1, const svec &v2);
  //! Template instantiation of elem_mult
  extern template bvec elem_mult(const bvec &v1, const bvec &v2);

  //------------- Elementwise Multiplication operator (three vectors) ----------

  //! Template instantiation of elem_mult
  extern template vec elem_mult(const vec &v1, const vec &v2, const vec &v3);
  //! Template instantiation of elem_mult
  extern template cvec elem_mult(const cvec &v1, const cvec &v2, const cvec &v3);
  //! Template instantiation of elem_mult
  extern template ivec elem_mult(const ivec &v1, const ivec &v2, const ivec &v3);
  //! Template instantiation of elem_mult
  extern template svec elem_mult(const svec &v1, const svec &v2, const svec &v3);
  //! Template instantiation of elem_mult
  extern template bvec elem_mult(const bvec &v1, const bvec &v2, const bvec &v3);

  //------------- Elementwise Multiplication operator (four vectors) ----------

  //! Template instantiation of elem_mult
  extern template vec elem_mult(const vec &v1, const vec &v2, const vec &v3, const vec &v4);
  //! Template instantiation of elem_mult
  extern template cvec elem_mult(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4);
  //! Template instantiation of elem_mult
  extern template ivec elem_mult(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4);
  //! Template instantiation of elem_mult
  extern template svec elem_mult(const svec &v1, const svec &v2, const svec &v3, const svec &v4);
  //! Template instantiation of elem_mult
  extern template bvec elem_mult(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4);

  //------------- Division operator ----------

  //! Template instantiation of operator/
  extern template vec operator/(const vec &v, double t);
  //! Template instantiation of operator/
  extern template cvec operator/(const cvec &v, std::complex<double> t);
  //! Template instantiation of operator/
  extern template ivec operator/(const ivec &v, int t);
  //! Template instantiation of operator/
  extern template svec operator/(const svec &v, short t);
  //! Template instantiation of operator/
  extern template bvec operator/(const bvec &v, bin t);

  //! Template instantiation of operator/
  extern template vec operator/(const double t, const vec &v);
  //! Template instantiation of operator/
  extern template cvec operator/(const std::complex<double> t, const cvec &v);
  //! Template instantiation of operator/
  extern template ivec operator/(const int t, const ivec &v);
  //! Template instantiation of operator/
  extern template svec operator/(const short t, const svec &v);
  //! Template instantiation of operator/
  extern template bvec operator/(const bin t, const bvec &v);

  //------------- Elementwise Division operator ----------

  //! Template instantiation of elem_div
  extern template vec elem_div(const vec &v1, const vec &v2);
  //! Template instantiation of elem_div
  extern template cvec elem_div(const cvec &v1, const cvec &v2);
  //! Template instantiation of elem_div
  extern template ivec elem_div(const ivec &v1, const ivec &v2);
  //! Template instantiation of elem_div
  extern template svec elem_div(const svec &v1, const svec &v2);
  //! Template instantiation of elem_div
  extern template bvec elem_div(const bvec &v1, const bvec &v2);

  //! Template instantiation of elem_div
  extern template vec elem_div(const double t, const vec &v);
  //! Template instantiation of elem_div
  extern template cvec elem_div(const std::complex<double> t, const cvec &v);
  //! Template instantiation of elem_div
  extern template ivec elem_div(const int t, const ivec &v);
  //! Template instantiation of elem_div
  extern template svec elem_div(const short t, const svec &v);
  //! Template instantiation of elem_div
  extern template bvec elem_div(const bin t, const bvec &v);

  //--------------------- concat operator -----------------

  //! Template instantiation of concat
  extern template vec concat(const vec &v, const double a);
  //! Template instantiation of concat
  extern template cvec concat(const cvec &v, const std::complex<double> a);
  //! Template instantiation of concat
  extern template ivec concat(const ivec &v, const int a);
  //! Template instantiation of concat
  extern template svec concat(const svec &v, const short a);
  //! Template instantiation of concat
  extern template bvec concat(const bvec &v, const bin a);

  //! Template instantiation of concat
  extern template vec concat(const double a, const vec &v);
  //! Template instantiation of concat
  extern template cvec concat(const std::complex<double> a, const cvec &v);
  //! Template instantiation of concat
  extern template ivec concat(const int a, const ivec &v);
  //! Template instantiation of concat
  extern template svec concat(const short a, const svec &v);
  //! Template instantiation of concat
  extern template bvec concat(const bin a, const bvec &v);

  //! Template instantiation of concat
  extern template vec concat(const vec &v1, const vec &v2);
  //! Template instantiation of concat
  extern template cvec concat(const cvec &v1, const cvec &v2);
  //! Template instantiation of concat
  extern template ivec concat(const ivec &v1, const ivec &v2);
  //! Template instantiation of concat
  extern template svec concat(const svec &v1, const svec &v2);
  //! Template instantiation of concat
  extern template bvec concat(const bvec &v1, const bvec &v2);

  //! Template instantiation of concat
  extern template vec concat(const vec &v1, const vec &v2, const vec &v3);
  //! Template instantiation of concat
  extern template cvec concat(const cvec &v1, const cvec &v2, const cvec &v3);
  //! Template instantiation of concat
  extern template ivec concat(const ivec &v1, const ivec &v2, const ivec &v3);
  //! Template instantiation of concat
  extern template svec concat(const svec &v1, const svec &v2, const svec &v3);
  //! Template instantiation of concat
  extern template bvec concat(const bvec &v1, const bvec &v2, const bvec &v3);

  //! Template instantiation of concat
  extern template vec concat(const vec &v1, const vec &v2, const vec &v3, const vec &v4);
  //! Template instantiation of concat
  extern template cvec concat(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4);
  //! Template instantiation of concat
  extern template ivec concat(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4);
  //! Template instantiation of concat
  extern template svec concat(const svec &v1, const svec &v2, const svec &v3, const svec &v4);
  //! Template instantiation of concat
  extern template bvec concat(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4);

  //! Template instantiation of concat
  extern template vec concat(const vec &v1, const vec &v2, const vec &v3, const vec &v4, const vec &v5);
  //! Template instantiation of concat
  extern template cvec concat(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4, const cvec &v5);
  //! Template instantiation of concat
  extern template ivec concat(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4, const ivec &v5);
  //! Template instantiation of concat
  extern template svec concat(const svec &v1, const svec &v2, const svec &v3, const svec &v4, const svec &v5);
  //! Template instantiation of concat
  extern template bvec concat(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4, const bvec &v5);

  // -------------- output stream --------------------

  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream& os, const vec &vect);
  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream& os, const cvec &vect);
  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream& os, const svec &vect);
  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream& os, const ivec &vect);
  //! Template instantiation of output stream
  extern template std::ostream &operator<<(std::ostream& os, const bvec &vect);

  // -------------- input stream --------------------

  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream& is, vec &vect);
  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream& is, cvec &vect);
  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream& is, svec &vect);
  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream& is, ivec &vect);
  //! Template instantiation of input stream
  extern template std::istream &operator>>(std::istream& is, bvec &vect);

#endif

} //namespace itpp

#endif // __vec_h
