
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
  \brief Definition of IT++ version function
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#ifndef __itpp_version_h
#define __itpp_version_h

#include <itpp/itconfig.h>
#include <string>


namespace itpp {

/*!
  \brief IT++ version function
  
  Returns the version number of IT++. E.g. "3.7.1". If your IT++ distribution
  is checked out from CVS "CVS version" is returned.
*/
  std::string itpp_version(void);

} //namespace itpp

#endif // __itpp_version_h
