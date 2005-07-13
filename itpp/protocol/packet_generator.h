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
  \brief Packet generator classes
  \author Anders Persson and Tony Ottosson

  $Revision$

  $Date$ 
*/

#ifndef __packet_generator_h
#define __packet_generator_h

#include "itpp/protocol/signals_slots.h"
#include "itpp/protocol/packet.h"
#include "itpp/base/random.h"

namespace itpp {

  /*!

  */
  class Packet_Generator {
  public:
    Packet_Generator(const int Packet_size = 150, const unsigned long int Max_packets = 0);
    virtual ~Packet_Generator();
    Signal<Packet*> output;
    Slot<Packet_Generator, bool> start;
    void set_parameters(const int Packet_size, const unsigned long int Max_packets);
    int get_packet_size();
    int get_max_packets();
  protected:
    virtual Ttype delta_t() = 0;
  private:
    Slot<Packet_Generator, Packet*> next;
    void handle_next(Packet*);
    void handle_start(const bool run);
    bool keep_running;
    unsigned long int id;
    int packet_size;
    unsigned long int max_packets;
  };


  /*!

  */
  class Poisson_Packet_Generator : public Packet_Generator {
  public:
    Poisson_Packet_Generator(const double Avg_bit_rate = 1.0, const int Packet_size = 150, const unsigned long int Max_packets = 0);
    virtual ~Poisson_Packet_Generator();
    void set_parameters(const double Avg_bit_rate, const int Packet_size, const unsigned long int Max_packets);
    double get_avg_bit_rate();
  protected:
    virtual Ttype delta_t();
    double avg_delta_t;
    double avg_bit_rate;
    Exponential_RNG ee;
  };


  /*!

  */
  class Constant_Rate_Packet_Generator : public Poisson_Packet_Generator {
  public:
    Constant_Rate_Packet_Generator(const double Avg_bit_rate = 1.0, const int Packet_size = 150, const unsigned long int Max_packets = 0);
    virtual ~Constant_Rate_Packet_Generator();
  protected:
    virtual Ttype delta_t();
  };

  /*!

  */
  class Burst_WWW_Packet_Generator : public Poisson_Packet_Generator {
  public:
    Burst_WWW_Packet_Generator(const double Avg_bit_rate = 1.0, const int Packet_size = 150, const int Max_packets = 0);
    virtual ~Burst_WWW_Packet_Generator();
  protected:
    virtual Ttype delta_t();
    int N, Navg;
    double Ti, Tr;
  };


  /*!

  */
  class Sink {
  public:
    Sink(const unsigned long int Max_packets = 1000);
    ~Sink();
    // -- Slots -- //
    Slot<Sink, Packet*> packet_input; 
  private:
    void handle_packet_input(Packet* packet);
    unsigned long int Ncp;
    unsigned long int Nbytes;
    unsigned long int max_packets;
    Ttype start_time;
  };

} // namespace itpp

#endif //__packet_generator_h

