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
  \brief Implementation of IT++ version function
  \author Tony Ottosson

  $Revision$

  $Date$
*/

#include "itpp/base/itpp_version.h"

namespace itpp { 

  std::string itpp_version(void)
  {
    std::string release_tag = "$Name$";
    std::string out;
    int pos;
  
    if (release_tag == "$Name$")
      return "CVS version";
    
    pos = release_tag.find_first_of("-");
    out = release_tag.substr(pos+1, release_tag.size());
  
    // Replace all remaining "-" with "."
    pos = out.find("-");
    while (pos != -1)
      {
	out.replace(pos, 1, ".");
	pos = out.find("-");
      }
  
    return out;
  }

} //namespace itpp
