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
  \brief Packet class
  \author Krister Norlund, Anders Persson and Tony Ottosson

  $Revision$

  $Date$ 
*/

#ifndef __packet_h
#define __packet_h

#include <stdlib.h>
#include <itpp/protocol/signals_slots.h>

namespace itpp {


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


} // namespace itpp

#endif //__packet_h

