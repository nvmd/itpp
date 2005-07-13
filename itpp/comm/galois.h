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
  \brief  Definitions of Galois Field algebra classes and functions
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __galois_h
#define __galois_h

#include "itpp/base/vec.h"
#include "itpp/base/specmat.h"
#include "itpp/base/array.h"
#include <iostream>
#include <cassert>

namespace itpp {

  /*! 
    \brief Galois Field GF(q).
    \author Tony Ottosson

    Galois field GF(q), where \a q = 2^m. Possible \a m values is \a m = 1,2,...,16.
    Elements are given as exponents of the primitive element \a alpha.
    Observe that the zeroth element are given as "-1". ( log(0)=-Inf ).
    <h3> The following primitve polynomials are used to construct the fields:</h3>
    <ul>
    <li> GF(4): 1+x+x^2 </li>
    <li> GF(8): 1+x+x^3 </li>
    <li> GF(16): 1+x+x^4 </li>
    <li> GF(32): 1+x^2+x^5 </li>
    <li> GF(64): 1+x^2+x^6 </li>
    <li> GF(128): 1+x^3+x^7 </li>
    <li> GF(256): 1+x^2+x^3+x^4+x^8 </li>
    <li> GF(512): 1+x^4+x^9 </li>
    <li> GF(1024): 1+x^3+x^10 </li>
    <li> GF(2^11): 1+x^2+x^11 </li>
    <li> GF(2^12): 1+x+x^4+x^12 </li>
    <li> GF(2^13): 1+x+x^3+x^4+x^13 </li>
    <li> GF(2^14): 1+x+x^3+x^5+x^14 </li>
    <li> GF(2^15): 1+x+x^15 </li>
    <li> GF(2^16): 1+x+x^3+x^12+x^16 </li>
    </ul>

    As indicated it is possible to use this class for binary elements, that is GF(2).
    However, this is less efficient
    in storage (each element take 5 bytes of memory) and in speed.
    If possible use the class BIN instead.
    Observe, also that the element "0" is called "-1" and "1" called "0".
  */
  class GF {
  public:
    //! Constructor
    GF() { m=0; }
    //! Constructor
    GF(int qvalue) { m=0; if (qvalue==0) // qvalue==0 gives the zeroth element
      value=-1; else set_size(qvalue); }
    //! Constructor
    GF(int qvalue, int inexp) { m=0; set(qvalue,inexp); }
    //! Copy constructor
    GF(const GF &ingf) { m=ingf.m; value=ingf.value; }

    //! GF(q) equals \a alpha ^ \a inexp
    void set(int qvalue, int inexp) {
      set_size(qvalue); it_assert0(inexp>=-1 && inexp<qvalue-1, "GF::set, out of range"); value=inexp; }
    /*! 
      \brief GF(q) equals the element that corresponds to the given vector space.

      The format is (...,c,b,a), where the element x is given as x=...+c*alpha^2+b*alpha+a.
    */
    void set(int qvalue, const bvec &vectorspace);
    //! set q=2^mvalue
    void set_size(int qvalue);
    //! Return q.
    int get_size() const { return ( (m != 0) ? q[m] : 0 ); }
    /*! 
      \brief Returns the vector space representation of GF(q).

      The format is (...,c,b,a), where the element x is given as x=...+c*alpha^2+b*alpha+a.
    */
    bvec get_vectorspace() const;
    //! Returns the alpha exponent
    int  get_value() const;
    //! Equality check
    int operator==(const GF &ingf) const;
    //! Not-equality check
    int operator!=(const GF &ingf) const;
  
    //! GF(q) equals ingf
    void operator=(const GF &ingf);
    //! GF(q) equals alpha^inexp
    void operator=(const int inexp);
    //! sum of two GF(q)
    void operator+=(const GF &ingf);
    //! sum of two GF(q)
    GF operator+(const GF &ingf) const;
    //! Difference of two GF(q), same as sum for q=2^m.
    void operator-=(const GF &ingf);
    //! Difference of two GF(q), same as sum for q=2^m.
    GF operator-(const GF &ingf) const;
    //! product of two GF(q)
    void operator*=(const GF &ingf);
    //! product of two GF(q)
    GF operator*(const GF &ingf) const;
    //! division of two GF(q)
    void operator/=(const GF &ingf);
    //! product of two GF(q)
    GF operator/(const GF &ingf) const;
    //! Output stream for GF(q)
    friend std::ostream &operator<<(std::ostream &os, const GF &ingf);
  protected:
  private:
    char m;
    int value;
    static Array<Array<int> > alphapow,logalpha;
    static ivec q;
  };

  class GFX;

  //! Multiplication of GF and GFX
  GFX  operator*(const GF &ingf, const GFX &ingfx);
  //! Multiplication of GFX and GF
  GFX  operator*( const GFX &ingfx, const GF &ingf);
  //! Division of GFX by GF
  GFX  operator/(const GFX &ingfx, const GF &ingf);

  /*! 
    \brief Polynomials over GF(q)[x], where q=2^m, m=1,...,16
  */
  class GFX {
  public:
    //! Constructor
    GFX();
    //! Constructor
    GFX(int qvalue);
    //! Constructor
    GFX(int qvalue, int indegree);
    //! Constructor
    GFX(int qvalue, const ivec &invalues);
    //! Constructor
    GFX(int qvalue, char *invalues);
    //! Constructor
    GFX(int qvalue, std::string invalues);
    //! Copy constructor
    GFX(const GFX &ingfx);
    //! Return q.
    int get_size() const;
    //! Return degree of GF(q)[x]
    int get_degree() const;
    /*! 
      \brief Resize the polynomial to the given indegree. If the new polynomial is bigger, then the new coefficients are set to zero.
    */
    void set_degree(int indegree);
    //! Return true degree of GF(q)[x]
    int get_true_degree() const;
    //! Set the GF(q)[x] polynomial
    void set(int qvalue, const char *invalues);
    //! Set the GF(q)[x] polynomial
    void set(int qvalue, const std::string invalues);
    //! Set the GF(q)[x] polynomial
    void set(int qvalue, const ivec &invalues);
    //! Set all coefficients to zero.
    void clear();
    //! Acces to individual element in the GF(q)[x] polynomial
    GF operator[](int index) const {
      it_assert0(index<=degree, "GFX::op[], out of range"); return coeffs(index); }
    //! Acces to individual element in the GF(q)[x] polynomial
    GF &operator[](int index) {
      it_assert0(index<=degree, "GFX::op[], out of range"); return coeffs(index); }
    //! Copy
    void operator=(const GFX &ingfx);
    //! sum of two GF(q)[x]
    void operator+=(const GFX &ingfx);
    //! sum of two GF(q)[x]
    GFX operator+(const GFX &ingfx) const;
    //! Difference of two GF(q), same as sum for q=2^m.
    void operator-=(const GFX &ingfx);
    //! Difference of two GF(q), same as sum for q=2^m.
    GFX operator-(const GFX &ingfx) const;
    //! product of two GF(q)[x]
    void operator*=(const GFX &ingfx);
    //! product of two GF(q)[x]
    GFX operator*(const GFX &ingfx) const;
    //! Evaluate polynom at alpha^inexp
    GF operator()(const GF &ingf);
    //! Multiply a GF element with a GF(q)[x]
    friend GFX  operator*(const GF &ingf, const GFX &ingfx);
    //! Multiply a GF(q)[x] with a GF element
    friend GFX  operator*( const GFX &ingfx, const GF &ingf);
    //! Divide a GF(q)[x] with a GF element
    friend GFX  operator/(const GFX &ingfx, const GF &ingf);

    //! Output stream
    friend std::ostream &operator<<(std::ostream &os, const GFX &ingfx);
  protected:
  private:
    int degree, q;
    Array<GF> coeffs;
  };

  //-------------- Help Functions ------------------
  /*! 
    \relates GFX
    \brief Int division of GF[q](x) polynomials: m(x) = c(x)/g(x).

    The reminder r(x) is not returned by this function.
  */
  GFX divgfx(const GFX &c, const GFX &g);

  /*! 
    \relates GFX
    \brief Function that performs int division of gf[q](x) polynomials (a(x)/g(x)) and returns the reminder.
  */
  GFX modgfx(const GFX &a, const GFX &b);


  // --------------- Inlines ------------------------
  // --------------- class GF -----------------------

  inline void GF::set(int qvalue, const bvec &vectorspace)
    {
      set_size(qvalue);
      it_assert0(vectorspace.length() == m, "GF::set, out of range");
      value=logalpha(m)(bin2dec(vectorspace));
    }

  inline bvec GF::get_vectorspace() const
    {
      bvec temp(m);
      if (value == -1)
	temp=dec2bin(m,0);
      else
	temp=dec2bin(m,alphapow(m)(value));
      return temp;
    }

  inline int  GF::get_value() const
    {
      return value;
    }

  inline int GF::operator==(const GF &ingf) const
    {
      if (value == -1 && ingf.value == -1)
	return true;
      if (m==ingf.m && value==ingf.value)
	return true;
      else
	return false;
    }

  inline int GF::operator!=(const GF &ingf) const
    {
      GF tmp(*this);
      return !(tmp==ingf);
    }

  inline void GF::operator=(const GF &ingf)
    {
      m=ingf.m;
      value=ingf.value;
    }

  inline void GF::operator=(const int inexp)
    {
      it_assert0(m>0 && inexp>=-1 && inexp<(q[m]-1), "GF::op=, out of range");
      value=inexp;
    }

  inline void GF::operator+=(const GF &ingf)
    {
      if (value == -1) {
	value=ingf.value;
	m=ingf.m;
      }
      else if (ingf.value != -1) {
	it_assert0(ingf.m == m, "GF::op+=, not same field");
	value=logalpha(m)(alphapow(m)(value)^alphapow(m)(ingf.value));
      }
    }

  inline GF GF::operator+(const GF &ingf) const
    {
      GF tmp(*this);
      tmp+=ingf;
      return tmp;
    }

  inline void GF::operator-=(const GF &ingf)
    {
      (*this)+=ingf;
    }

  inline GF GF::operator-(const GF &ingf) const
    {
      GF tmp(*this);
      tmp-=ingf;
      return tmp;
    }

  inline void GF::operator*=(const GF &ingf)
    {
      if (value == -1 || ingf.value == -1)
	value=-1;
      else {
	it_assert0(ingf.m == m, "GF::op+=, not same field");
	value=(value+ingf.value)%(q[m]-1);
      }
    }

  inline GF GF::operator*(const GF &ingf) const
    {
      GF tmp(*this);
      tmp*=ingf;
      return tmp;
    }

  inline void GF::operator/=(const GF &ingf)
    {
      assert(ingf.value !=-1); // no division by the zeroth element
      if (value == -1)
	value=-1;
      else {
	it_assert0(ingf.m == m, "GF::op+=, not same field");
	value=(value-ingf.value+q[m]-1)%(q[m]-1);
      }
    }

  inline GF GF::operator/(const GF &ingf) const
    {
      GF tmp(*this);
      tmp/=ingf;
      return tmp;
    }

  // ------------------ class GFX --------------------
  inline GFX::GFX()
    {
      degree=-1;
      q=0;
    }

  inline GFX::GFX(int qvalue)
    {
      it_assert0(qvalue>=0, "GFX::GFX, out of range");
      q=qvalue;
    }

  inline void GFX::set(int qvalue, const ivec &invalues)
    {
      it_assert0(qvalue>0, "GFX::set, out of range");
      degree=invalues.size()-1;
      coeffs.set_size(degree+1, false);
      for (int i=0;i<degree+1;i++)
	coeffs(i).set(qvalue,invalues(i));
      q=qvalue;
    }

  inline void GFX::set(int qvalue, const char *invalues)
    {
      set(qvalue,ivec(invalues));
    }

  inline void GFX::set(int qvalue, const std::string invalues)
    {
      set(qvalue,invalues.c_str());
    }

  inline GFX::GFX(int qvalue, int indegree)
    {
      it_assert0(qvalue>0 && indegree>=0, "GFX::GFX, out of range");
      q=qvalue;
      coeffs.set_size(indegree+1, false);
      degree=indegree;
      for (int i=0;i<degree+1;i++)
	coeffs(i).set(q,-1);
    }
  inline GFX::GFX(int qvalue, const ivec &invalues)
    {
      set(qvalue,invalues);
    }

  inline GFX::GFX(int qvalue, char *invalues)
    {
      set(qvalue,invalues);
    }

  inline GFX::GFX(int qvalue, std::string invalues)
    {
      set(qvalue,invalues.c_str());
    }

  inline GFX::GFX(const GFX &ingfx)
    {
      degree=ingfx.degree;
      coeffs=ingfx.coeffs;
      q=ingfx.q;
    }

  inline int GFX::get_size() const
    {
      return q;
    }

  inline int GFX::get_degree() const
    {
      return degree;
    }

  inline void GFX::set_degree(int indegree)
    {
      it_assert0(indegree>=-1, "GFX::set_degree, out of range");
      coeffs.set_size(indegree+1);
      degree=indegree;
    }

  inline int GFX::get_true_degree() const
    {
      int i=degree;
      while(coeffs(i).get_value()==-1) {
	i--;
	if (i==-1)
	  break;
      }
      return i;
    }

  inline void GFX::clear()
    {
      it_assert0(degree>=0 && q>0, "GFX::clear, not set");
      for(int i=0;i<degree+1;i++)
	coeffs(i).set(q,-1);
    }

  inline void GFX::operator=(const GFX &ingfx)
    {
      degree=ingfx.degree;
      coeffs=ingfx.coeffs;
      q=ingfx.q;
    }

  inline void GFX::operator+=(const GFX &ingfx)
    {
      it_assert0(q == ingfx.q, "GFX::op+=, not same field");
      if (ingfx.degree > degree) {
	coeffs.set_size(ingfx.degree+1, true);
	// set new coefficients to the zeroth element
	for (int j=degree+1; j<coeffs.size(); j++){ coeffs(j).set(q,-1); }
	degree=ingfx.degree;
      }
      for (int i=0;i<ingfx.degree+1;i++) { coeffs(i)+=ingfx.coeffs(i); }
    }

  inline GFX GFX::operator+(const GFX &ingfx) const
    {
      GFX tmp(*this);
      tmp+=ingfx;
      return tmp;
    }

  inline void GFX::operator-=(const GFX &ingfx)
    {
      (*this)+=ingfx;
    }

  inline GFX GFX::operator-(const GFX &ingfx) const
    {
      GFX tmp(*this);
      tmp-=ingfx;
      return tmp;
    }

  inline void GFX::operator*=(const GFX &ingfx)
    {
      it_assert0(q == ingfx.q, "GFX::op*=, Not same field");
      int i,j;
      Array<GF> tempcoeffs=coeffs;
      coeffs.set_size(degree+ingfx.degree+1, false);
      for (j=0; j<coeffs.size(); j++)
	coeffs(j).set(q,-1); // set coefficients to the zeroth element (log(0)=-Inf=-1)
      for (i=0;i<degree+1;i++)
	for (j=0;j<ingfx.degree+1;j++)
	  coeffs(i+j)+=tempcoeffs(i)*ingfx.coeffs(j);
      degree=coeffs.size()-1;
    }

  inline GFX GFX::operator*(const GFX &ingfx) const
    {
      GFX tmp(*this);
      tmp*=ingfx;
      return tmp;
    }

  inline GFX operator*(const GF &ingf, const GFX &ingfx)
    {
      it_assert0(ingf.get_size() == ingfx.q, "GFX::op*, Not same field");
      GFX temp(ingfx);
      for (int i=0;i<ingfx.degree+1;i++)
	temp.coeffs(i)*=ingf;
      return temp;
    }

  inline GFX  operator*( const GFX &ingfx, const GF &ingf)
    {
      return ingf*ingfx;
    }

  inline GFX  operator/(const GFX &ingfx, const GF &ingf)
    {
      it_assert0(ingf.get_size() == ingfx.q, "GFX::op/, Not same field");
      GFX temp(ingfx);
      for (int i=0;i<ingfx.degree+1;i++)
	temp.coeffs(i)/=ingf;
      return temp;
    }

  inline GF GFX::operator()(const GF &ingf)
    {
      it_assert0(q == ingf.get_size(), "GFX::op(), Not same field");
      GF temp(coeffs(0)), ingfpower(ingf);
      for (int i=1; i<degree+1; i++) {
	temp+=coeffs(i)*ingfpower;
	ingfpower*=ingf;
      }
      return temp;
    }

} //namespace itpp

#endif // __galois_h
