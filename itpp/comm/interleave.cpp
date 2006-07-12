/*!
 * \file 
 * \brief Implementation of interleaver classes
 * \author Pal Frenger
 *
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA
 *
 * -------------------------------------------------------------------------
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

} // namespace itpp
