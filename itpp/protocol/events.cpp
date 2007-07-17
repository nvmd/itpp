/*!
 * \file
 * \brief Implementation of an event-based simulation class
 * \author Anders Persson
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

#include <itpp/protocol/events.h>


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


} // namespace itpp
