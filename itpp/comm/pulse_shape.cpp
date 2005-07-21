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
  \brief Implementations of puls shaping classes.
  \author Tony Ottosson, Håkan Eriksson and Adam Piatyszek

  $Revision$

  $Date$
*/

#include <itpp/comm/pulse_shape.h>

namespace itpp {

  //---------------------------------------------------------------------------
  // Template instantiations
  //---------------------------------------------------------------------------

  template class Pulse_Shape<double, double, double>;
  template class Pulse_Shape<std::complex<double>, double, std::complex<double> >;
  template class Pulse_Shape<std::complex<double>, std::complex<double>, std::complex<double> >;

  template class Root_Raised_Cosine<double>;
  template class Root_Raised_Cosine<std::complex<double> >;

  template class Raised_Cosine<double>;
  template class Raised_Cosine<std::complex<double> >;
  
} //namespace itpp
