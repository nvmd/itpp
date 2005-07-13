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
  \brief Packet channel classes
  \author Anders Persson and Tony Ottosson

  $Revision$

  $Date$ 
*/

#ifndef __packet_channel_h
#define __packet_channel_h

#include "itpp/protocol/signals_slots.h"
#include "itpp/protocol/packet.h"
#include "itpp/base/vec.h"
#include "itpp/base/itassert.h"

namespace itpp {

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

} // namespace itpp

#endif //__packet_generator_h

