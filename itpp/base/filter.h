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
  \brief Definitions of Filter classes and functions
  \author Håkan Eriksson, Thomas Eriksson and Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __filter_h
#define __filter_h

#include <complex>
#include <cstdlib>
#include <iostream>
#include <algorithm>

#include "itpp/itconfig.h"
#include "itpp/base/vec.h"
#include "itpp/base/matfunc.h"
#include "itpp/base/specmat.h"
#include "itpp/base/itassert.h"

namespace itpp {

  /*!
    \addtogroup filters
  */

  /*!
    \brief Virtual Filter Base Class.
    \ingroup filters

    The class is templated as follows:
    <ul>
    <li> \c T1 is the type of the input samples</li>
    <li> \c T2 is the type of the filter coefficients</li>
    <li> \c T3 is the type of the output samples</li>
    </ul>
  */
  template <class T1, class T2, class T3>
    class Filter {
  public:
    //! Default constructor
    Filter();
    //! Filter a single sample.
    virtual T3 operator()(const T1 Sample) { return filter(Sample); }
    //! Filter a vector.
    virtual Vec<T3> operator()(const Vec<T1> &v);
  protected:
    /*!
      \brief Pure virtual filter function. This is where the real filtering is done. Implement this function to create a new filter.
    */
    virtual T3 filter(const T1 Sample)=0;
  };

  /*!
    \brief Moving Average Filter Base Class.
    \ingroup filters

    This class implements a moving average (MA) filter
    according to
    \f[
    y(n) = b(0)*x(n) + b(1)*x(n-1) + ... + b(N)*x(n-N)
    \f]
    where \a b is the filter coefficients, \a x is the input and \a y is the output.

    When filtering a vector, the length of the output vector equals the length of the input vector.
    Internal states are kept in a filter memory. The first time the filter is used the internal states
    have been set to zero.

    The class is templated as follows:
    <ul>
    <li> \c T1 is the type of the input samples</li>
    <li> \c T2 is the type of the filter coefficients</li>
    <li> \c T3 is the type of the output samples</li>
    </ul>
  */
  template <class T1, class T2, class T3>
    class MA_Filter : public Filter<T1,T2,T3> {
  public:
    //! Class default constructor
    explicit MA_Filter();
    //! Class constructor setting the coefficients in the filter
    explicit MA_Filter(const Vec<T2> &b);
    //! Class destructor
    virtual ~MA_Filter() { }
    //! Filter coefficient access function
    Vec<T2> get_coeffs() const { return coeffs; }
    //! Set the filter coefficients
    void set_coeffs(const Vec<T2> &b);
    //! Clears the filter memory
    void clear() { mem.clear(); }
    //! Get state of filter
    Vec<T3> get_state() const;
    //! Set state of filter
    void set_state(const Vec<T3> &state);

    //protected:
  private:
    virtual T3 filter(const T1 Sample);

    Vec<T3> mem;
    Vec<T2> coeffs;
    long inptr, NS;
    bool init;
  };

  /*!
    \brief Autoregressive (AR) Filter Base Class.
    \ingroup filters

    This class implements a autoregressive (AR) filter
    according to
    \f[
    a(0)*y(n) = x(n) - a(1)*y(n-1) - ... - a(N)*y(n-N)
    \f]
    where \a a is the filter coefficients, \a x is the input and \a y is the output.

    When filtering a vector, the length of the output vector equals the length of the input vector.
    Internal states are kept in a filter memory. The first time the filter is used the internal states
    have been set to zero.

    The class is templated as follows:
    <ul>
    <li> \c T1 is the type of the input samples</li>
    <li> \c T2 is the type of the filter coefficients</li>
    <li> \c T3 is the type of the output samples</li>
    </ul>
  */
  template <class T1, class T2, class T3>
    class AR_Filter : public Filter<T1,T2,T3> {
  public:
    //! Class constructor
    explicit AR_Filter();
    //! Class constructor setting the coefficients in the filter
    explicit AR_Filter(const Vec<T2> &a);
    //! Class destructor
    virtual ~AR_Filter() { }
    //! Filter coefficient access function
    Vec<T2> get_coeffs() const { return coeffs; }
    //! Set the filter coefficients (and order)
    void set_coeffs(const Vec<T2> &a);
    //! Clears the filter memory
    void clear() { mem.clear(); }
    //! Get state of filter
    Vec<T3> get_state() const;
    //! Set state of filter
    void set_state(const Vec<T3> &state);

    //protected:
  private:
    virtual T3 filter(const T1 Sample);

    Vec<T3> mem;
    Vec<T2> coeffs;
    long inptr, NS;
    bool init;
  };


  /*!
    \brief Autoregressive Moving Average (ARMA) Filter Base Class.
    \ingroup filters

    This class implements a autoregressive moving average (ARMA) filter
    according to
    \f[
    a(0)*y(n) = b(0)*x(n) + b(1)*x(n-1) + \ldots + b(N_b)*x(n-N_b) - a(1)*y(n-1) - \ldots - a(N_a)*y(n-N_a)
    \f]
    where \a a and \a b are the filter coefficients, \a x is the input and \a y is the output.

    When filtering a vector, the length of the output vector equals the length of the input vector.
    Internal states are kept in a filter memory. The first time the filter is used the internal states
    have been set to zero.

    The class is templated as follows:
    <ul>
    <li> \c T1 is the type of the input samples</li>
    <li> \c T2 is the type of the filter coefficients</li>
    <li> \c T3 is the type of the output samples</li>
    </ul>
  */
  template <class T1, class T2, class T3>
    class ARMA_Filter : public Filter<T1,T2,T3> {
  public:
    //! Class constructor
    explicit ARMA_Filter();
    //! Class constructor setting the coefficients in the filter
    explicit ARMA_Filter(const Vec<T2> &b, const Vec<T2> &a);
    //! Class destructor
    virtual ~ARMA_Filter() { }
    //! Filter \a a coefficient access function
    Vec<T2> get_coeffs_a() const { return acoeffs; }
    //! Filter \a b coefficient access function
    Vec<T2> get_coeffs_b() const { return bcoeffs; }
    //! Filter coefficient access function
    void get_coeffs(Vec<T2> &b, Vec<T2> &a) const { b = bcoeffs; a = acoeffs; }
    //! Set the filter coefficients (and order)
    void set_coeffs(const Vec<T2> &b, const Vec<T2> &a);
    //! Clears the filter memory
    void clear() { mem.clear(); }
    //! Get state of filter
    Vec<T3> get_state() const;
    //! Set state of filter
    void set_state(const Vec<T3> &state);

    //protected:
  private:
    virtual T3 filter(const T1 Sample);

    Vec<T3> mem;
    Vec<T2> acoeffs, bcoeffs;
    long inptr, NA, NB, NS;
    bool init;
  };



  /*!
    \brief ARMA filter function
    \ingroup filters

    These functions implements a autoregressive moving average (ARMA) filter
    according to
    \f[
    a(0)*y(n) = b(0)*x(n) + b(1)*x(n-1) + \ldots + b(N_b)*x(n-N_b) - a(1)*y(n-1) - \ldots - a(N_a)*y(n-N_a)
    \f]
    where \a a and \a b are the filter coefficients, \a x is the input and \a y is the output.
    Setting a=1 gives a MA filter and b=1 gives a AR filter.

    The length of the output vector equals the length of the input vector. The state vectors 
    \a state_in and \a state_out is of length \f$max(N_a, n_b) - 1\f$.
    
    If no start state \a state_in is given it is set to zero.
  */
  //@{
  vec filter(const vec &b, const vec &a, const vec &input);
  cvec filter(const vec &b, const vec &a, const cvec &input);
  cvec filter(const cvec &b, const cvec &a, const cvec &input);
  cvec filter(const cvec &b, const cvec &a, const vec &input);

  vec filter(const vec &b, const int one, const vec &input);
  cvec filter(const vec &b, const int one, const cvec &input);
  cvec filter(const cvec &b, const int one, const cvec &input);
  cvec filter(const cvec &b, const int one, const vec &input);
  
  vec filter(const int one, const vec &a, const vec &input);
  cvec filter(const int one, const vec &a, const cvec &input);
  cvec filter(const int one, const cvec &a, const cvec &input);
  cvec filter(const int one, const cvec &a, const vec &input);


  vec filter(const vec &b, const vec &a, const vec &input, const vec &state_in, vec &state_out);
  cvec filter(const vec &b, const vec &a, const cvec &input, const cvec &state_in, cvec &state_out);
  cvec filter(const cvec &b, const cvec &a, const cvec &input, const cvec &state_in, cvec &state_out);
  cvec filter(const cvec &b, const cvec &a, const vec &input, const cvec &state_in, cvec &state_out);

  vec filter(const vec &b, const int one, const vec &input, const vec &state_in, vec &state_out);
  cvec filter(const vec &b, const int one, const cvec &input, const cvec &state_in, cvec &state_out);
  cvec filter(const cvec &b, const int one, const cvec &input, const cvec &state_in, cvec &state_out);
  cvec filter(const cvec &b, const int one, const vec &input, const cvec &state_in, cvec &state_out);
  
  vec filter(const int one, const vec &a, const vec &input, const vec &state_in, vec &state_out);
  cvec filter(const int one, const vec &a, const cvec &input, const cvec &state_in, cvec &state_out);
  cvec filter(const int one, const cvec &a, const cvec &input, const cvec &state_in, cvec &state_out);
  cvec filter(const int one, const cvec &a, const vec &input, const cvec &state_in, cvec &state_out);
  //@}


  /*!
    \brief Design a Nth order FIR filter with cut-off frequency \c cutoff using the window method.
    \ingroup filters
  */
  vec fir1(long N, double cutoff);

  //-------------------------------------------------------------------------------------------------
  // Implementation of templated functions starts here
  //-------------------------------------------------------------------------------------------------

  //---------------------- class Filter ----------------------------

  template <class T1, class T2, class T3>
    Filter<T1,T2,T3>::Filter()
    {
    }

  template <class T1, class T2, class T3>
    Vec<T3> Filter<T1,T2,T3>::operator()(const Vec<T1> &x)
    {
      Vec<T3> y(x.length());

      for (long i=0;i<x.length();i++) {
	y[i]=filter(x[i]);
      }
      return y;
    }

  //-------------------------- class MA_Filter ---------------------------------

  template <class T1, class T2,class T3>
    MA_Filter<T1,T2,T3>::MA_Filter() : Filter<T1,T2,T3>()
    {
      inptr=0;
      init = false;
    }
    
    template <class T1, class T2,class T3>
      MA_Filter<T1,T2,T3>::MA_Filter(const Vec<T2> &b) : Filter<T1,T2,T3>()
      {
	set_coeffs(b);
      }


      template <class T1, class T2, class T3>
	void MA_Filter<T1,T2,T3>::set_coeffs(const Vec<T2> &b)
	{
	  it_assert(b.size() > 0, "MA_Filter: size of filter is 0!");

	  coeffs = b;
	  NS = coeffs.size();

	  mem.set_size(NS, false);
	  mem.clear();
	  inptr=0;
	  init = true;
	}

      template <class T1, class T2, class T3>
	Vec<T3> MA_Filter<T1,T2,T3>::get_state() const
	{
	  it_assert(init == true, "MA_Filter: filter coefficients are not set!");

	  int offset = inptr;
	  Vec<T3> state(NS);

	  for(int n = 0; n<NS; n++){
	    state(n) = mem(offset);
	    offset = (offset+1)%NS;
	  }

	  return state;
	}

      template <class T1, class T2, class T3>
	void MA_Filter<T1,T2,T3>::set_state(const Vec<T3> &state)
	{
	  it_assert(init == true, "MA_Filter: filter coefficients are not set!");
	  it_assert(state.size() == NS, "MA_Filter: Invalid state vector!");

	  mem = state;
	  inptr = 0;
	}

      template <class T1, class T2, class T3>
	T3 MA_Filter<T1,T2,T3>::filter(const T1 Sample)
	{
	  it_assert(init == true, "MA_Filter: Filter coefficients are not set!");
	  T3  s=0;
	  long i,L;

	  mem._elem(inptr)=Sample;
	  L=mem.length()-inptr;

	  for (i=0;i<L;i++) {
	    s+=coeffs._elem(i)*mem._elem(inptr+i);
	  }
	  for (i=0;i<inptr;i++) {
	    s+=coeffs._elem(L+i)*mem._elem(i);
	  }
	  inptr=inptr-1;if (inptr<0) inptr+=mem.length();
	  return s;
	}

      //---------------------- class AR_Filter ----------------------------------

      template <class T1, class T2, class T3>
	AR_Filter<T1,T2,T3>::AR_Filter() : Filter<T1,T2,T3>()
	{
	  inptr = 0;
	  init = false;
	}

	template <class T1, class T2, class T3>
	  AR_Filter<T1,T2,T3>::AR_Filter(const Vec<T2> &a) : Filter<T1,T2,T3>()
	  {
	    set_coeffs(a);
	  }

	  template <class T1, class T2, class T3>
	    void AR_Filter<T1,T2,T3>::set_coeffs(const Vec<T2> &a)
	    {
	      it_assert(a.size() > 0, "AR_Filter: size of filter is 0!");
	      it_assert(a(0) != T2(0), "AR_Filter: a(0) cannot be 0!");

	      coeffs = a.right(length(a)-1)/a(0);
	      NS = coeffs.size();

	      mem.set_size(NS, false);
	      mem.clear();
	      inptr = 0;
	      init = true;
	    }


	  template <class T1, class T2, class T3>
	    Vec<T3> AR_Filter<T1,T2,T3>::get_state() const
	    {
	      it_assert(init == true, "AR_Filter: filter coefficients are not set!");

	      int offset = inptr;
	      Vec<T3> state(NS);

	      for(int n = 0; n<NS; n++){
		state(n) = mem(offset);
		offset = (offset+1)%NS;
	      }

	      return state;
	    }

	  template <class T1, class T2, class T3>
	    void AR_Filter<T1,T2,T3>::set_state(const Vec<T3> &state)
	    {
	      it_assert(init == true, "AR_Filter: filter coefficients are not set!");
	      it_assert(state.size() == NS, "AR_Filter: Invalid state vector!");

	      mem = state;
	      inptr = 0;
	    }

	  template <class T1, class T2, class T3>
	    T3 AR_Filter<T1,T2,T3>::filter(const T1 Sample)
	    {
	      it_assert(init == true, "AR_Filter: Filter coefficients are not set!");
	      T3  s=0;
	      long i,L;

	      L=length(mem)-inptr;
	      for (i=0;i<L;i++) {
		s+=coeffs._elem(i)*mem._elem(inptr+i);
	      }
	      for (i=0;i<inptr;i++) {
		s+=coeffs._elem(L+i)*mem._elem(i);
	      }
	      inptr=inptr-1;if (inptr<0) inptr+=length(mem);

	      s=Sample-s;
	      mem._elem(inptr)=s;

	      return s;
	    }


	  //---------------------- class ARMA_Filter ----------------------------------
	  template <class T1, class T2, class T3>
	    ARMA_Filter<T1,T2,T3>::ARMA_Filter() : Filter<T1,T2,T3>()
	    {
	      inptr = 0;
	      init = false;
	    }

	    template <class T1, class T2, class T3>
	      ARMA_Filter<T1,T2,T3>::ARMA_Filter(const Vec<T2> &b, const Vec<T2> &a) : Filter<T1,T2,T3>()
	      {
		set_coeffs(b, a);
	      }

	      template <class T1, class T2, class T3>
		void ARMA_Filter<T1,T2,T3>::set_coeffs(const Vec<T2> &b, const Vec<T2> &a)
		{
		  it_assert(a.size() > 0 && b.size() > 0, "ARMA_Filter: size of filter is 0!");
		  it_assert(a(0) != T2(0), "ARMA_Filter: a(0) cannot be 0!");

		  acoeffs = a/a(0);
		  bcoeffs = b/a(0);

		  inptr = 0;
		  NA = a.size();
		  NB = b.size();
		  NS = std::max(NB, NA) - 1;
		  mem.set_size(NS, false);
		  mem.clear();
		  init = true;
		}

	      template <class T1, class T2, class T3>
		Vec<T3> ARMA_Filter<T1,T2,T3>::get_state() const
		{
		  it_assert(init == true, "ARMA_Filter: filter coefficients are not set!");

		  int offset = inptr;
		  Vec<T3> state(NS);

		  for(int n = 0; n<NS; n++){
		    state(n) = mem(offset);
		    offset = (offset+1)%NS;
		  }
      
		  return state;
		}

	      template <class T1, class T2, class T3>
		void ARMA_Filter<T1,T2,T3>::set_state(const Vec<T3> &state)
		{
		  it_assert(init == true, "ARMA_Filter: filter coefficients are not set!");
		  it_assert(state.size() == NS, "ARMA_Filter: Invalid state vector!");

		  mem = state;
		  inptr = 0;
		}

	      template <class T1, class T2, class T3>
		T3 ARMA_Filter<T1,T2,T3>::filter(const T1 Sample)
		{
		  it_assert(init == true, "ARMA_Filter: Filter coefficients are not set!");
		  const int N0 = NS-1;
		  int ia, ib;
		  T3 z = Sample;
		  T3 s;
      
		  for(ia=0; ia<NA-1; ia++) // All AR-coeff except a(0).
		    z -= mem((ia+inptr)%NS) * acoeffs(ia+1);

		  s = z*bcoeffs(0);

		  for(ib=0; ib<NB-1; ib++) // All MA-coeff except b(0).
		    s += mem((ib+inptr)%NS) * bcoeffs(ib+1);

		  inptr = (inptr+N0)%NS; // Clock the state shift-register once.
		  mem(inptr) = z; // Store in the internal state.

		  return s;
		}



	      //-----------------------------------------------------------------------
	      //  class Filter
	      //-----------------------------------------------------------------------
#ifndef _MSC_VER


	      //-----------------------------------------------------------------------
	      //  class MA_Filter
	      //-----------------------------------------------------------------------

	      //! Template instatiation of MA_Filter
	      extern template class MA_Filter<double,double,double>;
	      //! Template instatiation of MA_Filter
	      extern template class MA_Filter<double,std::complex<double>,std::complex<double> >;
	      //! Template instatiation of MA_Filter
	      extern template class MA_Filter<std::complex<double>,double,std::complex<double> >;
	      //! Template instatiation of MA_Filter
	      extern template class MA_Filter<std::complex<double>,std::complex<double>,std::complex<double> >;

	      //-----------------------------------------------------------------------
	      //  class AR_Filter
	      //-----------------------------------------------------------------------

	      //! Template instatiation of AR_Filter
	      extern template class AR_Filter<double,double,double>;
	      //! Template instatiation of AR_Filter
	      extern template class AR_Filter<double,std::complex<double>,std::complex<double> >;
	      //! Template instatiation of AR_Filter
	      extern template class AR_Filter<std::complex<double>,double,std::complex<double> >;
	      //! Template instatiation of AR_Filter
	      extern template class AR_Filter<std::complex<double>,std::complex<double>,std::complex<double> >;

	      //-----------------------------------------------------------------------
	      //  class ARMA_Filter
	      //-----------------------------------------------------------------------

	      //! Template instatiation of AR_Filter
	      extern template class ARMA_Filter<double,double,double>;
	      //! Template instatiation of AR_Filter
	      extern template class ARMA_Filter<double,std::complex<double>,std::complex<double> >;
	      //! Template instatiation of AR_Filter
	      extern template class ARMA_Filter<std::complex<double>,double,std::complex<double> >;
	      //! Template instatiation of AR_Filter
	      extern template class ARMA_Filter<std::complex<double>,std::complex<double>,std::complex<double> >;


#endif // MSC_VER

} //namespace itpp

#endif // __filter_h
