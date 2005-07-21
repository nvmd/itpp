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
  \brief Class factory for fixed-point data types Fix and CFix
  \author Johan Bergman
  
  $Revision$
  
  $Date$
*/

#include <cstring>
#include <itpp/fixedpoint/fix_factory.h>
#include <itpp/fixedpoint/fix.h>
#include <itpp/fixedpoint/cfix.h>

namespace itpp {

  void Fix_Factory::create(Fix* &ptr, const int n) const
  {
    ptr = new Fix[n];
    if (ptr) {
      // Set fixed-point restrictions
      for (int i=0; i<n; i++) {
        // Note: explicit destructor call intentionally omitted (to save time)
        new (ptr + i) Fix(0.0, 0, wordlen, emode, omode, qmode, stat_ptr);
      }
    }
  }

  void Fix_Factory::create(CFix* &ptr, const int n) const
  {
    ptr = new CFix[n];
    if (ptr) {
      // Set fixed-point restrictions
      for (int i=0; i<n; i++) {
        // Note: explicit destructor call intentionally omitted (to save time)
        new (ptr + i) CFix(0.0, 0, wordlen, emode, omode, qmode, stat_ptr);
      }
    }
  }

  template<>
  void create_elements<Fix>(Fix* &ptr, const int n, const Factory &f)
  {
    if (const Fix_Factory *fix_factory_ptr = dynamic_cast<const Fix_Factory*>(&f)) {
      // Yes, f seems to be a Fix_Factory. Now call the Fix_Factory::create method
      fix_factory_ptr->create(ptr, n);
    }
    else {
      // No, f does not seem to be a Fix_Factory. As a fallback solution,
      // assume that f is DEFAULT_FACTORY and use the default constructor
      ptr = new Fix[n];
    }
  }

  template<>
  void create_elements<CFix>(CFix* &ptr, const int n, const Factory &f)
  {
    if (const Fix_Factory *fix_factory_ptr = dynamic_cast<const Fix_Factory*>(&f)) {
      // Yes, f seems to be a Fix_Factory. Now call the Fix_Factory::create method
      fix_factory_ptr->create(ptr, n);
    }
    else {
      // No, f does not seem to be a Fix_Factory. As a fallback solution,
      // assume that f is DEFAULT_FACTORY and use the default constructor
      ptr = new CFix[n];
    }
  }

} //namespace itpp
