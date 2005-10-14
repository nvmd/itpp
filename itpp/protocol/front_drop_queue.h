/*!
 * \file 
 * \brief Definitions of a Front Drop Queue class
 * \author Anders Persson and Tony Ottosson
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

#ifndef FRONT_DROP_QUEUE_H
#define FRONT_DROP_QUEUE_H

#include <queue>
#include <itpp/protocol/signals_slots.h>
#include <itpp/protocol/events.h>
#include <itpp/protocol/packet.h>

namespace itpp {

#define DEFAULT_MAX_BYTES_IN_QUEUE 24000

class Front_Drop_Queue : public virtual std::queue<Packet*> {
  public:
   Front_Drop_Queue(const int max_bytes = DEFAULT_MAX_BYTES_IN_QUEUE)  {
      max_bytes_in_queue = max_bytes;
      bytes_in_queue = 0;
      debug=false;
   }

   // TODO destructor
//  ~FrontDropQueue() { }

  void set_debug(const bool enable_debug = true) { 
     debug = enable_debug;
  }

  void push(Packet *packet);
  void pop();

  void set_max_byte_size(int max_bytes) { max_bytes_in_queue = max_bytes; }
  int max_byte_size() { return max_bytes_in_queue; }
  int byte_size() { return bytes_in_queue; }

 private:
  int max_bytes_in_queue;
  int bytes_in_queue;
  int debug;
};



} // namespace itpp

#endif // #ifndef FRONT_DROP_QUEUE_H

