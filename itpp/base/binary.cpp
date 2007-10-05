/*!
 * \file
 * \brief Binary class implemenations
 * \author Tony Ottosson
 *
 * -------------------------------------------------------------------------
 *
 * IT++ - C++ library of mathematical, signal processing, speech processing,
 *        and communications classes and functions
 *
 * Copyright (C) 1995-2007  (see AUTHORS file for a list of contributors)
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

#include <itpp/base/binary.h>
#include <iostream>


namespace itpp {

  std::ostream &operator<<(std::ostream &output, const bin &inbin)
  {
    output << static_cast<int>(inbin);
    return output;
  }

  std::istream &operator>>(std::istream &input, bin &outbin)
  {
    int tmp;
    input >> tmp;
    it_assert((tmp == 0) || (tmp == 1),
              "bin::operator>>(): input value must be 0 or 1");
    outbin = tmp;
    return input;
  }

} // namespace itpp
