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
  \brief Binary class implemenations
  \author Tony Ottosson

  $Revision$

  $Date$
*/
 
#include "itpp/base/binary.h"

using std::ostream;
using std::istream;

namespace itpp { 

  ostream &operator<<(ostream &output, const bin &inbin)
  {
    output << static_cast<int>(inbin);
    return output;
  }

  istream &operator>>(istream &input, bin &outbin)
  {
    int tmp;
    input >> tmp;
    outbin = tmp;
    return input;
  }

} //namespace itpp
