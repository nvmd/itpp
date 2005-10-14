/*!
 * \file
 * \brief Implementation of IT++ version function
 * \author Tony Ottosson
 * 
 * $Date$
 * $Revision$
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2005  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/itpp_version.h>

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
    while (pos != -1) {
			out.replace(pos, 1, ".");
			pos = out.find("-");
		}
		
    return out;
  }

} //namespace itpp
