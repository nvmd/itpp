/*!
 * \file 
 * \brief Definition of a Packet class
 * \author Krister Norlund, Anders Persson and Tony Ottosson
 *
 * $Date$
 * $Revision$
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

#ifndef PACKET_H
#define PACKET_H

#include <itpp/protocol/signals_slots.h>


namespace itpp {

  //! \addtogroup protocol
  //@{

  /*! \brief Packet

  */
  class Packet {
  public:
    //!
    Packet(const int packet_size=0) { set_bit_size(packet_size); }
    //!
    virtual ~Packet() {}

    //! set size of packet in bits
    void set_bit_size(int packet_size) { it_assert(packet_size >= 0, "Packet size must be positive"); size_bits = packet_size; }

    //! get size of packet in bits
    int bit_size() { return size_bits; }

  private:
    int size_bits; // size of packet in bits
  };


  /*!

  */
  class L3_Packet_Info{
  public:
    L3_Packet_Info(Packet *packet) { timestamp = 0; pkt_pointer = packet; }

    ~L3_Packet_Info() {}

    Ttype timestamp;

    Packet *pkt_pointer;
  };


  /*!

  */
  class Link_Packet : public Packet {
  public:
    Link_Packet(const int Seq_no, const unsigned long int Link_packet_id, L3_Packet_Info *Cp) { seq_no = Seq_no; link_packet_id = Link_packet_id; l3_pkt_info_p = Cp; }

    ~Link_Packet() {}

    unsigned long int link_packet_id;
    int seq_no;
    L3_Packet_Info *l3_pkt_info_p;
  };

  /*!

  */
  class ACK : public Packet {
  public:
    ACK(const int Seq_no=-1, const int Id=0) { seq_no = Seq_no; id = Id; }

    ~ACK() {}

    int id;
    int seq_no;
  };

  //@}

} // namespace itpp

#endif // #ifndef PACKET_H

