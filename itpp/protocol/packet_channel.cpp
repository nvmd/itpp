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

#include <itpp/protocol/packet_channel.h>
#include <itpp/base/random.h>
#include <itpp/base/sort.h>
#include <itpp/base/stat.h>

namespace itpp {

Packet_Channel::Packet_Channel() {
  parameters_ok = false;
  keep_running = false;
}

Packet_Channel::Packet_Channel(const double Pr, const Ttype Delay, const double Block_rate, const int Max_slots) {
  set_parameters(Pr, Delay, Block_rate, Max_slots);
}


Packet_Channel::~Packet_Channel(){}

void Packet_Channel::set_parameters(const double Pr,  const Ttype Delay,  const double Block_rate, const int Max_slots) {
  assert(Delay >= 0);
  assert(Pr>=0.0 && Pr<=1.0);
  assert(Block_rate > 0);
  assert(Max_slots >= 0);
  delay = Delay;
  pr = Pr;
  block_time = 1.0/Block_rate;
  max_slots = Max_slots;
  input.forward(this, &Packet_Channel::handle_input);
  nof_inputs.forward(this, &Packet_Channel::handle_nof_inputs);
  start.forward(this, &Packet_Channel::handle_start);
  keep_running = false;
  explicit_errors = false;
  K = 0;
  k = 0;
  parameters_ok = true;
}

void Packet_Channel::handle_input(Link_Packet* M) {
  assert(parameters_ok);
  assert(M!=NULL);
  if(explicit_errors){
    if(k<L){
      lose = lost(k)==K;
      if(lose)
	k++;
    }
    K++;
  }
  else
    lose = randu() < pr;
  if(lose){
    delete M;  
  }
  else
    output(M, delay);  
  lose = false;
}

void Packet_Channel::block_rate_loop() {
  assert(parameters_ok);
  get_nof_inputs(NULL);
  if(keep_running)
    Event_Queue::add(new Event<Packet_Channel>(this, &Packet_Channel::block_rate_loop, block_time));
}

void Packet_Channel::handle_start(const bool run) {
  assert(parameters_ok);
  if(run&&!keep_running)// Channel is in 'stop' state. Start it and keep running.
    Event_Queue::add(new Event<Packet_Channel>(this, &Packet_Channel::block_rate_loop, block_time));  
  keep_running = run;
}

void Packet_Channel::handle_nof_inputs(const int Nof_ready_messages) {
  assert(Nof_ready_messages>=0);
  int L = 0;
  if(max_slots>0)
    L = std::min(Nof_ready_messages, round_i(randu()*max_slots));
  else
    L = std::min(Nof_ready_messages, 1);
  if(L>0)
    input_request(L);
}

void Packet_Channel::set_errors(const ivec &Lost) {
  L = Lost.length();
  if(L>0){
    assert(min(Lost)>=0);
    lost = Lost;
    sort(lost);
    explicit_errors = true;
  }
}


  // ----------------------------- Ack_Channel --------------------------------


ACK_Channel::ACK_Channel() {
  parameters_ok = false;
}

ACK_Channel::ACK_Channel(const double Pr, const Ttype Delay) {
  set_parameters(Pr, Delay);
}


ACK_Channel::~ACK_Channel(){}

void ACK_Channel::set_parameters(const double Pr, const Ttype Delay) {
  assert(Delay >= 0);
  assert(Pr>=0.0 && Pr<=1.0);
  delay = Delay;
  pr = Pr;
  input.forward(this, &ACK_Channel::handle_input);
  explicit_errors = false;
  K = 0;
  k = 0;
  parameters_ok = true;
}

void ACK_Channel::handle_input(ACK* M) {
  assert(parameters_ok);
  assert(M!=NULL);
  if(explicit_errors){
    if(k<L){
      lose = lost(k)==K;
      if(lose)
	k++;
    }
    K++;
  }
  else
    lose = randu() < pr;
  if(lose)
    delete M;
  else
    output(M, delay);  
  lose = false;
}

void ACK_Channel::set_errors(const ivec& Lost) {
  L = Lost.length();
  if(L>0){
    assert(min(Lost)>=0);
    lost = Lost;
    sort(lost);
    explicit_errors = true;
  }
}



} //namespace itpp
