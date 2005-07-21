/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 1995-2001 by Tony Ottosson, Thomas Eriksson, Pål Frenger,   *
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
  \brief Definitions of deterministic sources
  \author Tobias Ringström and Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __source_h
#define __source_h

#include <cstdlib>

#include <itpp/base/binary.h>
#include <itpp/base/vec.h>
#include <itpp/base/mat.h>
#include <itpp/base/scalfunc.h>

namespace itpp {

  /*! 
    \addtogroup detsource
  */

  /*! 
    \brief A sine-wave source
    \ingroup detsource
  */
  class Sine_Source {
  public:
    //! Constructor. Set frequency, mean, amplitude, and start phase
    Sine_Source(double freq, double mean=0.0, double ampl=1.0, double inphase=0.0);
    //! Return a single sample
    double operator()() { return sample(); }
    //! Get a sample vector.
    vec operator()(int n);
    //! Get a sample matrix.
    mat operator()(int h, int w);
  protected:
  private:
    double sample();
    double m, A, theta, dtheta;
  };

  /*! 
    \brief A square-wave source
    \ingroup detsource
  */
  class Square_Source {
  public: 
    //! Constructor. Set frequency, mean, amplitude, and start phase
    Square_Source(double freq, double mean=0.0, double ampl=1.0, double inphase=0.0);
    //! Return a single sample
    double operator()() { return sample(); }
    //! Get a sample vector.
    vec operator()(int n);
    //! Get a sample matrix.
    mat operator()(int h, int w);
  protected:
  private:
    double sample();
    double m, A, theta, dtheta;
  };

  /*! 
    \brief A triangle-wave source
    \ingroup detsource
  */
  class Triangle_Source {
  public:
    //! Constructor. Set frequency, mean, amplitude and start phase
    Triangle_Source(double freq, double mean=0.0, double ampl=1.0, double inphase=0.0);
    //! Return a single sample
    double operator()() { return sample(); }
    //! Get a sample vector.
    vec operator()(int n);
    //! Get a sample matrix.
    mat operator()(int h, int w);
  protected:
  private:
    double sample();
    double m, A, theta, dtheta;
  };

  /*! 
    \brief A sawtooth-wave source
    \ingroup detsource
  */
  class Sawtooth_Source {
  public:
    //! Constructor. Set frequency, mean, amplitude, and start phase
    Sawtooth_Source(double freq, double mean=0.0, double ampl=1.0, double inphase=0.0);
    //! Return a single sample
    double operator()() { return sample(); }
    //! Get a sample vector.
    vec operator()(int n);
    //! Get a sample matrix.
    mat operator()(int h, int w);
  protected:
  private:
    double sample();
    double m, A, theta, dtheta;
  };

  /*! 
    \brief An impulse source
    \ingroup detsource
  */
  class Impulse_Source {
  public:
    //! Constructor. Set frequency, amplitude and start phase
    Impulse_Source(double freq, double ampl=1.0, double inphase=0.0);
    //! Return a single sample
    double operator()() { return sample(); }
    //! Get a sample vector.
    vec operator()(int n);
    //! Get a sample matrix.
    mat operator()(int h, int w);
  protected:
  private:
    double sample();
    double A, pos, dtheta;
  };

  /*! 
    \brief An pattern source
    \ingroup detsource
  */
  class Pattern_Source {
  public:
    //! Constructor. Set pattern and start position
    Pattern_Source(const vec &pattern, int start_pos=0);
    //! Destructor
    virtual ~Pattern_Source() { }
    //! Return a single sample
    double operator()() { return sample(); }
    //! Get a sample vector.
    vec operator()(int n);
    //! Get a sample matrix.
    mat operator()(int h, int w);
  protected:
  private:
    double sample();
    int pos;
    vec pat;
    double mean, var;
  };

} //namespace itpp

#endif // __source_h
