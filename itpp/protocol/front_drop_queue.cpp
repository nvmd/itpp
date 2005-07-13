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
  \brief Front Drop Queue class
  \author Anders Persson and Tony Ottosson

  $Revision$

  $Date$ 
*/

#include "itpp/protocol/front_drop_queue.h"


namespace itpp {


void Front_Drop_Queue::push(Packet *packet)
{
   if (debug) {
      std::cout << "Front_Drop_Queue::push_packet"
//                << " byte_size=" << packet->bit_size()/8
                << " ptr=" << packet
                << " time=" << Event_Queue::now() << std::endl;
   }
   
   Packet *hol_packet;
   while ((!std::queue<Packet*>::empty()) && 
          ((8*bytes_in_queue + packet->bit_size()) >  8*max_bytes_in_queue)) {
      hol_packet = std::queue<Packet*>::front();
      Front_Drop_Queue::pop();
      delete hol_packet;

//      TTCPPacket *tcp_packet = (TTCPPacket *) hol_packet;
//      delete tcp_packet;

      if (debug) {
         std::cout << "Link_With_Input_Q::received_packet, "
                   << "Packet Dropped, buffer overflow."
                   << std::endl;
      }
   }     

   bytes_in_queue += packet->bit_size()/8;
   std::queue<Packet*>::push(packet);

}

void Front_Drop_Queue::pop()
{
   Packet *hol_packet;
   hol_packet = std::queue<Packet*>::front();
   bytes_in_queue -= (hol_packet->bit_size()/8);
   if (debug) {
      std::cout << "Front_Drop_Queue::pop_packet"
                << " ptr=" << hol_packet
                << " time=" << Event_Queue::now() << std::endl;
   }
   std::queue<Packet*>::pop();

}

} //namespace itpp
