/*!
 * \file
 * \brief Definition of a Packet channel class
 * \author Anders Persson and Tony Ottosson
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

#ifndef PACKET_CHANNEL_H
#define PACKET_CHANNEL_H

#include <itpp/protocol/packet.h>
#include <itpp/base/vec.h>


namespace itpp {

  //! \addtogroup protocol
  //@{

  /*!

  */
  class Packet_Channel {
  public:
    Packet_Channel();
    Packet_Channel(const double Pr, const Ttype Delay, const double Block_rate, const int Max_slots = 0);

    ~Packet_Channel();

    // -- Slots -- //
    Slot<Packet_Channel, bool> start;
    Slot<Packet_Channel, Link_Packet*> input;
    Slot<Packet_Channel, int> nof_inputs;

    // -- Signals -- //
    Signal<Link_Packet*> output;
    Signal<int> input_request;
    Signal<void*> get_nof_inputs;

    void set_parameters(const double Pr, const Ttype Delay, const double Block_rate, const int Max_slots);

    void set_errors(const ivec &Lost);

  private:
    void block_rate_loop();
    void handle_input(Link_Packet* M);
    void handle_start(const bool start);
    void handle_nof_inputs(const int N);

    bool keep_running;
    bool parameters_ok;
    bool explicit_errors;
    bool lose;
    double pr;
    Ttype delay;
    double block_time;
    int max_slots;
    ivec lost;
    int k,K,L;
  };


  /*!

  */
  class ACK_Channel {
  public:
    ACK_Channel();

    ACK_Channel(const double Pr, const Ttype Delay);

    ~ACK_Channel();

    // -- Slots -- //
    Slot<ACK_Channel, ACK*> input;

    // -- Signals -- //
    Signal<ACK*> output;

    void set_parameters(const double Pr, const Ttype Delay);
    void set_errors(const ivec& Lost);

  private:
    void handle_input(ACK* M);

    bool parameters_ok;
    bool explicit_errors;
    bool lose;
    double pr;
    Ttype delay;
    ivec lost;
    int k, K, L;
  };

  //@}

} // namespace itpp

#endif // #ifndef PACKET_CHANNEL_H

