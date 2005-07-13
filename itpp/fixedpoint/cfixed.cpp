/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 2005 by Johan Bergman.                                      *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Templated complex fixed-point data type CFixed
  \author Johan Bergman
  
  $Revision$
  
  $Date$
*/

#include "itpp/fixedpoint/cfixed.h"

namespace itpp {

  // Template instantiations
  template class CFixed<64, TC, WRAP>;

} //namespace itpp
