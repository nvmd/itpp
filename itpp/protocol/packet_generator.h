/*!
 * \file 
 * \brief Definition of a Packet generator class
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
 * Copyright (C) 1995-2006  (see AUTHORS file for a list of contributors)
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

#ifndef PACKET_GENERATOR_H
#define PACKET_GENERATOR_H

#include <itpp/protocol/packet.h>
#include <itpp/base/random.h>


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

#endif // #ifndef PACKET_GENERATOR_H

