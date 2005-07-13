/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2004 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Implementation of deterministic sources.
  \author Tobias Ringström and Tony Ottosson

  $Revision$

  $Date$
*/

#include <ctime>
//#ifdef _MSC_VER
//#  include <process.h>
//#  define getpid _getpid
//#else
//#  include <unistd.h>
//#endif

#include "itpp/base/binary.h"
#include "itpp/base/matfunc.h"
#include "itpp/base/scalfunc.h"
#include "itpp/base/source.h"

namespace itpp { 

#ifndef DOXYGEN_SHOULD_SKIP_THIS

#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif
#ifndef M_PI_4
#define M_PI_4 0.78539816339744830962
#endif

#endif //DOXYGEN_SHOULD_SKIP_THIS

  ///////////////////////////////////////////////
  // Sine_Source
  ///////////////////////////////////////////////

  Sine_Source::Sine_Source(double freq, double mean, double ampl, double inphase)
  {
    A = ampl;
    m = mean;
    theta = inphase;
    dtheta = 2.0 * pi * freq;
  }

  double Sine_Source::sample()
  {
    double samp = m + A * sin(theta);
    
    theta += dtheta;
    if (theta >= 2.0 * pi)
      theta -= 2.0 * pi;

    return samp;
  }

  vec Sine_Source::operator()(int n)
  {
    vec v(n);

    for (int i=0; i<n; i++)
      v(i) = sample();

    return v;
  }

  mat Sine_Source::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // Square_Source
  ///////////////////////////////////////////////

  Square_Source::Square_Source(double freq, double mean, double ampl, double inphase)
  {
    A = ampl;
    m = mean;
    theta = inphase / (2.0 * pi);
    dtheta = freq;
  }

  double Square_Source::sample()
  {
    double samp = theta < 0.5 ? 1.0 : -1.0;
    
    theta += dtheta;
    if (theta >= 1.0)
      theta -= 1.0;

    return samp;
  }

  vec Square_Source::operator()(int n)
  {
    vec v(n);

    for (int i=0; i<n; i++)
      v(i) = sample();

    return v;
  }

  mat Square_Source::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // Triangle_Source
  ///////////////////////////////////////////////

  Triangle_Source::Triangle_Source(double freq, double mean, double ampl, double inphase)
  {
    A = ampl;
    m = mean;
    theta = inphase / (2.0 * pi);
    dtheta = freq;
  }

  double Triangle_Source::sample()
  {
    double samp = m + 4.0 * A * (theta < 0.25 ? theta : 0.5 - theta);
    
    theta += dtheta;
    if (theta >= 0.75)
      theta -= 1.0;

    return samp;
  }

  vec Triangle_Source::operator()(int n)
  {
    vec v(n);

    for (int i=0; i<n; i++)
      v(i) = sample();

    return v;
  }

  mat Triangle_Source::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // Sawtooth_Source
  ///////////////////////////////////////////////

  Sawtooth_Source::Sawtooth_Source(double freq, double mean, double ampl, double inphase)
  {
    A = ampl;
    m = mean;
    theta = inphase / (2.0 * pi);
    dtheta = freq;
  }

  double Sawtooth_Source::sample()
  {
    double samp = 2.0 * A * theta;
    
    theta += dtheta;
    if (theta >= 0.5)
      theta -= 1.0;

    return samp;
  }

  vec Sawtooth_Source::operator()(int n)
  {
    vec v(n);

    for (int i=0; i<n; i++)
      v(i) = sample();

    return v;
  }

  mat Sawtooth_Source::operator()(int h, int w)
  {
    mat mm(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	mm(i,j) = sample();

    return mm;
  }

  ///////////////////////////////////////////////
  // Impulse_Source
  ///////////////////////////////////////////////

  Impulse_Source::Impulse_Source(double freq, double ampl, double inphase)
  {
    A = ampl;
    pos = inphase / (2.0 * pi);
    dtheta = freq;
  }

  double Impulse_Source::sample()
  {
    double samp;
    
    if (pos >= 1.0) {
      samp = A;
      pos -= 1.0;
    }
    else {
      samp = 0.0;
      pos += dtheta;
    }

    return samp;
  }

  vec Impulse_Source::operator()(int n)
  {
    vec v(n);

    for (int i=0; i<n; i++)
      v(i) = sample();

    return v;
  }

  mat Impulse_Source::operator()(int h, int w)
  {
    mat m(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	m(i,j) = sample();

    return m;
  }

  ///////////////////////////////////////////////
  // Pattern_Source
  ///////////////////////////////////////////////

  Pattern_Source::Pattern_Source(const vec &pattern, int start_pos)
  {
    pat = pattern;
    pos = start_pos;

    // Calculate the mean and variance.  Note that the variance shall
    // be normalied by N and not N-1 in this case
    mean = var = 0.0;
    for (int i=pat.size()-1; i>=0; i--) {
      mean += pat(i);
      var += pat(i) * pat(i);
    }
    mean /= pat.size();
    var /= pat.size();
    var -= mean*mean;
  }

  double Pattern_Source::sample()
  {
    double samp = pat(pos);
    
    if (pos >= pat.size()-1)
      pos = 0;
    else
      pos++;

    return samp;
  }

  vec Pattern_Source::operator()(int n)
  {
    vec v(n);

    for (int i=0; i<n; i++)
      v(i) = sample();

    return v;
  }

  mat Pattern_Source::operator()(int h, int w)
  {
    mat m(h,w);
    int i,j;

    for (i=0; i<h; i++)
      for (j=0; j<w; j++)
	m(i,j) = sample();

    return m;
  }

  // -------------------------------------------------------------------------------

} //namespace itpp
