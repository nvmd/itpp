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
  \brief Implementation of Filter classes and functions.
  \author Håkan Eriksson, Thomas Eriksson and Tony Ottosson

  $Revision$

  $Date$
*/

#include "itpp/base/filter.h"

namespace itpp {



  vec filter(const vec &b, const vec &a, const vec &input)
  { 
    ARMA_Filter<double, double, double> f(b, a);
    return f(input);
  }

  cvec filter(const vec &b, const vec &a, const cvec &input)
  {
    ARMA_Filter<std::complex<double>,double,std::complex<double> > f(b, a);
    return f(input);
  }

  cvec filter(const cvec &b, const cvec &a, const cvec &input)
  {
    ARMA_Filter<std::complex<double>,std::complex<double>,std::complex<double> > f(b, a);
    return f(input);
  }

  cvec filter(const cvec &b, const cvec &a, const vec &input)
  {
    ARMA_Filter<double,std::complex<double>,std::complex<double> > f(b, a);
    return f(input);
  }


  vec filter(const vec &b, const int one, const vec &input)
  {
    it_assert(one == 1, "filter(): in a MA filter a=1");
    MA_Filter<double, double, double> f(b);
    return f(input);
  }

  cvec filter(const vec &b, const int one, const cvec &input)
  {
    it_assert(one == 1, "filter(): in a MA filter a=1");
    MA_Filter<std::complex<double>,double,std::complex<double> > f(b);
    return f(input);
  }

  cvec filter(const cvec &b, const int one, const cvec &input)
  {
    it_assert(one == 1, "filter(): in a MA filter a=1");
    MA_Filter<std::complex<double>,std::complex<double>,std::complex<double> > f(b);
    return f(input); }

  cvec filter(const cvec &b, const int one, const vec &input)
  {
    it_assert(one == 1, "filter(): in a MA filter a=1");
    MA_Filter<double,std::complex<double>,std::complex<double> > f(b);
    return f(input);
  }


  vec filter(const int one, const vec &a, const vec &input)
  {
    it_assert(one == 1, "filter(): in a AR filter b=1");
    AR_Filter<double, double, double> f(a);
    return f(input);
  }

  cvec filter(const int one, const vec &a, const cvec &input)
  {
    it_assert(one == 1, "filter(): in a AR filter b=1");
    AR_Filter<std::complex<double>,double,std::complex<double> > f(a);
    return f(input);
  }

  cvec filter(const int one, const cvec &a, const cvec &input)
  {
    it_assert(one == 1, "filter(): in a AR filter b=1");
    AR_Filter<std::complex<double>,std::complex<double>,std::complex<double> > f(a);
    return f(input);
  }

  cvec filter(const int one, const cvec &a, const vec &input)
  {
    it_assert(one == 1, "filter(): in a AR filter b=1");
    AR_Filter<double,std::complex<double>,std::complex<double> > f(a);
    return f(input);
  }





  vec filter(const vec &b, const vec &a, const vec &input, const vec &state_in, vec &state_out)
  { 
    ARMA_Filter<double, double, double> f(b, a);
    f.set_state(state_in);
    vec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const vec &b, const vec &a, const cvec &input, const cvec &state_in, cvec &state_out)
  {
    ARMA_Filter<std::complex<double>,double,std::complex<double> > f(b, a);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const cvec &b, const cvec &a, const cvec &input, const cvec &state_in, cvec &state_out)
  {
    ARMA_Filter<std::complex<double>,std::complex<double>,std::complex<double> > f(b, a);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const cvec &b, const cvec &a, const vec &input, const cvec &state_in, cvec &state_out)
  {
    ARMA_Filter<double,std::complex<double>,std::complex<double> > f(b, a);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }


  vec filter(const vec &b, const int one, const vec &input, const vec &state_in, vec &state_out)
  {
    it_assert(one == 1, "filter(): in a MA filter a=1");
    MA_Filter<double, double, double> f(b);
    f.set_state(state_in);
    vec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const vec &b, const int one, const cvec &input, const cvec &state_in, cvec &state_out)
  {
    it_assert(one == 1, "filter(): in a MA filter a=1");
    MA_Filter<std::complex<double>,double,std::complex<double> > f(b);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const cvec &b, const int one, const cvec &input, const cvec &state_in, cvec &state_out)
  {
    it_assert(one == 1, "filter(): in a MA filter a=1");
    MA_Filter<std::complex<double>,std::complex<double>,std::complex<double> > f(b);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const cvec &b, const int one, const vec &input, const cvec &state_in, cvec &state_out)
  {
    it_assert(one == 1, "filter(): in a MA filter a=1");
    MA_Filter<double,std::complex<double>,std::complex<double> > f(b);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }


  vec filter(const int one, const vec &a, const vec &input, const vec &state_in, vec &state_out)
  {
    it_assert(one == 1, "filter(): in a AR filter b=1");
    AR_Filter<double, double, double> f(a);
    f.set_state(state_in);
    vec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const int one, const vec &a, const cvec &input, const cvec &state_in, cvec &state_out)
  {
    it_assert(one == 1, "filter(): in a AR filter b=1");
    AR_Filter<std::complex<double>,double,std::complex<double> > f(a);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const int one, const cvec &a, const cvec &input, const cvec &state_in, cvec &state_out)
  {
    it_assert(one == 1, "filter(): in a AR filter b=1");
    AR_Filter<std::complex<double>,std::complex<double>,std::complex<double> > f(a);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }

  cvec filter(const int one, const cvec &a, const vec &input, const cvec &state_in, cvec &state_out)
  {
    it_assert(one == 1, "filter(): in a AR filter b=1");
    AR_Filter<double,std::complex<double>,std::complex<double> > f(a);
    f.set_state(state_in);
    cvec output = f(input);
    state_out = f.get_state();
    return output;
  }

  vec fir1(long N, double cutoff)
  {
    vec a(N+1),h=hamming(N+1);

    for (long i=0;i<length(a);i++) {
      a[i]=h[i]*sinc(cutoff*(i-N/2.0));
    }
    a/=sum(a);
    return a;
  }

  //-----------------------------------------------------------------------
  //  class Filter
  //-----------------------------------------------------------------------

  //template class Filter<double,double,double>;
  //template class Filter<double,std::complex<double>,std::complex<double> >;
  //template class Filter<std::complex<double>,double,std::complex<double> >;
  //template class Filter<std::complex<double>,std::complex<double>,std::complex<double> >;

  //-----------------------------------------------------------------------
  //  class MA_Filter
  //-----------------------------------------------------------------------

  template class MA_Filter<double,double,double>;
  template class MA_Filter<double,std::complex<double>,std::complex<double> >;
  template class MA_Filter<std::complex<double>,double,std::complex<double> >;
  template class MA_Filter<std::complex<double>,std::complex<double>,std::complex<double> >;

  //-----------------------------------------------------------------------
  //  class AR_Filter
  //-----------------------------------------------------------------------

  template class AR_Filter<double,double,double>;
  template class AR_Filter<double,std::complex<double>,std::complex<double> >;
  template class AR_Filter<std::complex<double>,double,std::complex<double> >;
  template class AR_Filter<std::complex<double>,std::complex<double>,std::complex<double> >;

  //-----------------------------------------------------------------------
  //  class ARMA_Filter
  //-----------------------------------------------------------------------

  template class ARMA_Filter<double,double,double>;
  template class ARMA_Filter<double,std::complex<double>,std::complex<double> >;
  template class ARMA_Filter<std::complex<double>,double,std::complex<double> >;
  template class ARMA_Filter<std::complex<double>,std::complex<double>,std::complex<double> >;



} //namespace itpp
