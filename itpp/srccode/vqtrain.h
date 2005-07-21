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
  \brief Vector quantizer trainer
  \author Thomas Eriksson

  $Revision$

  $Date$
*/

#ifndef __vqtrain_h
#define __vqtrain_h

#include <itpp/base/mat.h>
#include <itpp/base/array.h>

namespace itpp {


  //!
  double kmeansiter(Array<vec> &DB, mat &codebook);
  //!
  mat kmeans(Array<vec> &DB, int SIZE, int NOITER=9999, bool VERBOSE=true);
  //!
  mat lbg(Array<vec> &DB, int SIZE, int NOITER=9999, bool VERBOSE=true);


  /*! 
    \ingroup sourcecoding
    \brief Function for vector quantization training

    The following code illustrates how the VQ can be trained.

    \code
    VQ Quantizer;
    mat A;
    Array<vec> database;

    // read vectors into database somehow
    ...

    // train a vq
    A = vqtrain(database, 1024, 1000000);
    Quantizer.set_codebook(A);
    \endcode
  */
  mat vqtrain(Array<vec> &DB, int SIZE, int NOITER, double STARTSTEP=0.2, bool VERBOSE=true);
  //!
  vec sqtrain(const vec &inDB, int SIZE);

  //!
  ivec bitalloc(const vec& variances, int nobits); 
} // namespace itpp

#endif
