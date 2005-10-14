/*!
 * \file 
 * \brief Implementation of a class factory for fixed-point data types Fix 
 * and CFix
 * \author Johan Bergman
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

} // namespace itpp
