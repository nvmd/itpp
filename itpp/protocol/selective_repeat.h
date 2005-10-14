/*!
 * \file 
 * \brief Definitions of Selective Repeat ARQ classes
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

#ifndef SELECTIVE_REPEAT_H
#define SELECTIVE_REPEAT_H

#include <stdlib.h>
#include <itpp/protocol/signals_slots.h>
#include <itpp/protocol/packet.h>
#include <itpp/protocol/front_drop_queue.h>
#include <itpp/base/array.h>
#include <itpp/base/vec.h>

namespace itpp {

  /*!

  */
class Selective_Repeat_ARQ_Sender {
public:
  Selective_Repeat_ARQ_Sender();

  Selective_Repeat_ARQ_Sender(const int Seq_no_size, const int Buffer_size_factor, const int Link_packet_size, const Ttype Time_out);

  ~Selective_Repeat_ARQ_Sender();

  void set_parameters(const int Seq_no_size,        // # bits in sequence no. 
		      const int Buffer_size_factor, // Link-packet buffer size = 2^(Seq_no_size)*Buffer_size_factor.
		      const int Link_packet_size,    // Size of the link packets in bytes.
		      const Ttype Time_out);        // Idle time before retransmission.

  // -- Slots -- //
  Slot<Selective_Repeat_ARQ_Sender, Packet*> packet_input;            // Receives incoming packets.
  Slot<Selective_Repeat_ARQ_Sender, Array<Packet*> > ack_input;                  // Receives incoming ack/nacks.
  Slot<Selective_Repeat_ARQ_Sender, void*> query_nof_ready_packets;   // Receives incoming query for number of packets ready to transmit.
  Slot<Selective_Repeat_ARQ_Sender, int> packet_output_request;       // Receives incoming packet output requests.

  // -- Signals -- //
  Signal<Array<Packet*> > packet_output;                // Delivers transmitted packets.
  Signal<int> nof_ready_packets;                // Delivers no ready packets.
  Signal<int> buffer_overflow;                  // Signals buffer overflows.
  
  int buffer_size();
  int link_packets_buffered();
  int nof_ready_link_packets();
  int link_packets_queued_waiting_for_transmission();
  Ttype link_packets_max_queuing_time();
  void get_link_packets(const int K, Array<Packet*> &pa);

private:
  void handle_ack_input(Array<Packet*> packet_array); // Take care of incomming ack/nacks.
  void handle_packet_input(Packet *P);          // Take care of incomming packets.
  void handle_packet_output_request(int K);     // Take care of incomming packet requests.
  void handle_query_nof_ready_packets(void*);   // Take care of incomming query for number of packets ready to transmit.
  void retransmit(int Sequence_number);    // Take care of incomming query for number of packets ready to transmit.
  void remove(const int Sequence_number);
  void push_packet_on_tx_buffer(Packet *packet);
  int buffered_non_outstanding();
  int free_sequence_numbers();
  int sequence_number_2_buffer_index(const int Sequence_number);
  void schedule_output(const int Buffer_index, const int Sequence_number, const bool Retransmission);
  void cancel_output(const int Sequence_number);
  void fill_output();
  int feasable_blocks();
  bool parameters_ok;
  Front_Drop_Queue ip_pkt_queue;
  Array<Link_Packet*> input_buffer;
  int input_buffer_size;
  int input_next;
  int input_free_space;  
  int seq_no_size;
  int seq_no;
  int seq_no_max;
  int tx_next;
  int tx_last;
  int outstanding;
  int id;
  Ttype time_out;
  Array<ATimer<Selective_Repeat_ARQ_Sender, int> > timer;
  ivec output_indexes;
  ivec retransmission_indexes;
  int rd_pos;
  int rt_pos;
  int scheduled_total;
  int scheduled_retransmissions;
  int no_retransmit;
  int link_packet_size;
};


/*!

*/
class Selective_Repeat_ARQ_Receiver {
public:
  Selective_Repeat_ARQ_Receiver();

  Selective_Repeat_ARQ_Receiver(const int Seq_no_size);

  ~Selective_Repeat_ARQ_Receiver();

  // -- Slots -- //
  Slot<Selective_Repeat_ARQ_Receiver, Array<Packet*> > packet_input;    // Receives incoming packets.

  // -- Signals -- //
  Signal<Array<Packet*> > ack_output;       // Delivers ack.  
  Signal<Packet*> packet_output; // Delivers received packets.  

  void set_parameters(const int Seq_no_size); // # bits in sequence no. 

private:
  bool greater_modulo_L(const int a, const int b);
  void handle_packet_input(Array<Packet*>); // Take care of incomming packets.
  int seq_no_size;
  int seq_no_max;
  Array<Link_Packet*> rx_buffer;
  int Rnext;
  int id;
  bool parameters_ok;
};


} // namespace itpp

#endif // #ifndef SELECTIVE_REPEAT_H

