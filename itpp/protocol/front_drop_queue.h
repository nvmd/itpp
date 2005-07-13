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
  \brief Front Drop Queue
  \author Anders Persson and Tony Ottosson

  $Revision$

  $Date$ 
*/


#ifndef __front_drop_queue_h
#define __front_drop_queue_h


#include <queue>
#include "itpp/protocol/signals_slots.h"
#include "itpp/protocol/events.h"
#include "itpp/protocol/packet.h"


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

#endif //__front_drop_queue_h

