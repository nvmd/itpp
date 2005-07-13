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
  \brief Event-based simulation
  \author Anders Persson

  $Revision$

  $Date$ 
*/

#include "itpp/protocol/events.h"


namespace itpp {

  Ttype Event_Queue::t = 0;

  unsigned long long int Base_Event::global_id = 0;

  std::priority_queue<Base_Event*, 
		      std::deque<Base_Event*, std::allocator<Base_Event*> >,
		      Compare_Base_Event_Times> Event_Queue::event_queue;

  bool Event_Queue::keep_running = false;

  void Event_Queue::add(Base_Event *e)
  {
    e->expire_t = t + e->delta_t;
    event_queue.push(e);  
  }

  void Event_Queue::_run()
  {
    while(!event_queue.empty() && keep_running) {
      Base_Event* e = event_queue.top(); // Next event to process.
      event_queue.pop(); // Remove event.

      if(e->active) { // Only process active events.
	t = e->expire_t; // Update current time.
	e->exec(); // Execute the event.		
      }

      delete e; // This event is history!
    }

  }

  void Event_Queue::start()
  {
    keep_running = true;
    _run();
  }

  void Event_Queue::stop()
  {
    keep_running = false;  
  }

  void Event_Queue::clear()
  {
    stop();
    Base_Event* e; 

    while(!event_queue.empty()) {
      e = event_queue.top();
      delete e;
      event_queue.pop();
    }

    t = 0;
  }

  // void Event_Queue::cancel_all(BaseSignal *s){
  
  // }


} //namespace itpp
