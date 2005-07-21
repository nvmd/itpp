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
  \brief Implementation of Interleaver classes
  \author Pål Frenger

  $Revision$

  $Date$
*/

#include <itpp/comm/interleave.h>

namespace itpp {

  // ------------------------ Instantiations --------------------------------

  template class Block_Interleaver<double>;
  template class Block_Interleaver<short>;
  template class Block_Interleaver<int>;
  template class Block_Interleaver<std::complex<double> >;
  template class Block_Interleaver<bin>;

  template class Cross_Interleaver<double>;
  template class Cross_Interleaver<short>;
  template class Cross_Interleaver<int>;
  template class Cross_Interleaver<std::complex<double> >;
  template class Cross_Interleaver<bin>;

  template class Sequence_Interleaver<double>;
  template class Sequence_Interleaver<short>;
  template class Sequence_Interleaver<int>;
  template class Sequence_Interleaver<std::complex<double> >;
  template class Sequence_Interleaver<bin>;

} //namespace itpp
