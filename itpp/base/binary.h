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
  \brief Binary class definitions 
  \author Tony Ottosson and Sumit Bhardwaj

  $Revision$

  $Date$
*/

#ifndef __binary_h
#define __binary_h

#include <iostream>
#include <cassert>

#include <itpp/base/itassert.h>

namespace itpp {

  /*! 
    \brief Binary arithmetic (boolean) class
    \author Tony Ottosson

    This class creates a binary aritmetic class, following the ordinary
    rules for binary (GF(2)) fields.

    Examples:
    \code
    bin a;         // Creation of variable
    bin a = 0;     // Creating a variable and assigning it value 0
    bin b = 1;     // Creating a variable and assigning it value 1
    bin c = a + b; // XOR operation
    c = !a;        // NOT
    c = a * b;     // AND
    c = a / b;     // OR
    \endcode
  */
  class bin {
  public:
    //! Default constructor
    bin() { b=0; }
    
    //! Set the binary object equal to \c value. Either "0" or "1".
    bin(const short value) {
      it_assert0(value==0 || value==1, "bin(value): value must be 0 or 1");
      b = static_cast<char>(value);
    }

    //! Copy constructor
    bin(const bin &inbin) { b=inbin.b; }

    //! Assign a value
    void operator=(const short &value) {
      it_assert0(value==0 || value==1, "bin(value): value must be 0 or 1");
      b = static_cast<char>(value);
    }

    //! Assign a value
    void operator=(const bin &inbin)   { b=inbin.b; }

    //! OR
    void operator/=(const bin &inbin) { b=b|inbin.b; }

    //! OR
    void operator|=(const bin &inbin) { b=b|inbin.b; }
    //! OR
    bin operator/(const bin &inbin) const { return bin(b|inbin.b); }
    //! OR
    bin operator|(const bin &inbin) const { return bin(b|inbin.b); }

    //! XOR
    void operator+=(const bin &inbin) { b=b^inbin.b; }
    //! XOR
    void operator^=(const bin &inbin) { b=b^inbin.b; }
    //! XOR
    bin operator+(const bin &inbin) const { return bin(b^inbin.b); }
    //! XOR
    bin operator^(const bin &inbin) const { return bin(b^inbin.b); }
    //! XOR
    void operator-=(const bin &inbin) { b=b^inbin.b; }
    //! XOR
    bin operator-(const bin &inbin) const {return bin(b^inbin.b); }
    //! Dummy definition to be able to use vec<bin>
    bin operator-() const { return bin(b); }
	
    //! AND
    void operator*=(const bin &inbin) { b=b&inbin.b; }
    //! AND
    void operator&=(const bin &inbin) { b=b&inbin.b; }
    //! AND
    bin operator*(const bin &inbin) const { return bin(b&inbin.b); }
    //! AND
    bin operator&(const bin &inbin) const { return bin(b&inbin.b); }
	
    //! NOT
    bin operator!(void) const { return bin(b^1); }
    //! NOT
    bin operator~(void) const { return bin(b^1); }

    //! Check if equal
    bool operator==(const bin &inbin) const { return b == inbin.b; }
    //! Check if equal
    bool operator==(const int &inbin) const { return b == inbin; }

    //! Check if not equal
    bool operator!=(const bin &inbin) const { return b != inbin.b; }
    //! Check if not equal
    bool operator!=(const int &inbin) const { return b != inbin; }

    //! Less than (interpret the binary values {0,1} as integers)
    bool operator<(const bin &inbin) const  { return b < inbin.b; }
    //! Less than equal (interpret the binary values {0,1} as integers)
    bool operator<=(const bin &inbin) const { return b <= inbin.b; }

    //! Greater than (interpret the binary values {0,1} as integers)
    bool operator>(const bin &inbin) const  { return b > inbin.b; }
    //! Greater than equal (interpret the binary values {0,1} as integers)
    bool operator>=(const bin &inbin) const { return b >= inbin.b; }
	
    //! Convert \c bin to \c short
    operator short() const  { return static_cast<short>(b); }
    //! Convert \c bin to \c int
    operator int() const    { return static_cast<int>(b); }
    //! Convert \c bin to \c bool
    operator bool() const   { return b!=0; }
    //! Convert \c bin to \c float
    operator float() const  { return static_cast<float>(b); }
    //! Convert \c bin to \c double
    operator double() const { return static_cast<double>(b); }

    //! Output the binary value of the object
    char value() const { return b; }

  protected:

  private:
    char b;
  };



  /*! 
    \relates bin
    \brief Output stream of bin
  */
  std::ostream &operator<<(std::ostream &output, const bin &inbin);

  /*! 
    \relates bin
    \brief Input stream of bin
  */
  std::istream &operator>>(std::istream &input, bin &outbin);


  /*! 
    \relates bin
    \brief absolute value of bin
  */
  inline bin abs(const bin &inbin) { return inbin; }

} //namespace itpp

#endif // __binary_h


