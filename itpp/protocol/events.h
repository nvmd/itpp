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


#ifndef __events_h
#define __events_h

#include <list>
#include <queue>
#include <deque>
#include <string>
#include <iostream>

#include <itpp/base/itassert.h>

namespace itpp {
  //typedef long double Ttype; // 128-bit floating point time.
  typedef double Ttype; // 64-bit floating point time.
  //typedef long unsigned int Ttype; // 64-bit unsigned integer time.

  class Event_Queue;
  class Base_Event;
  class Base_Signal;

  /*!

  */
  class Base_Event {
  public:
    friend class Base_Signal;

    friend class Event_Queue;

    friend struct Compare_Base_Event_Times;
  
    //! Schedule an event at time \c delta_time from now
    Base_Event(const Ttype delta_time) {  // The event will occur in 'delta_time' time units from now! 	  
      it_assert(delta_time>=0, "Only causal simulations are possible");
      active = true;
      delta_t = delta_time;
      expire_t = 0; // Will be set correctly upon calling Event_Queue::add().
      id = global_id++;
    }

    //!
    virtual ~Base_Event(){}

    //! Cancel an event
    void cancel(){ active = false; } 

  protected:
    virtual void exec(void) = 0;
    Ttype delta_t;
    Ttype expire_t;
    bool active;
    unsigned long long int id;
    static unsigned long long int global_id;
  };

  //!
  struct Compare_Base_Event_Times {
    //!
    bool operator()(Base_Event *event1, Base_Event *event2) {
      if(event1->expire_t == event2->expire_t) // Equal expire times.
	return (event1->id > event2->id); // Base comparison on the event id.   
      else
	return (event1->expire_t > event2->expire_t); // Different expire times. Regular comparison.
    }
  };

  /*!

  */
  class Event_Queue {
  public:
    friend class Base_Signal;

    //!
    Event_Queue(){}
    //!
    ~Event_Queue(){}
    
    //! Add event to Queue
    static void add(Base_Event *e);  
    //! Return current time
    static Ttype now(){return t;}  
    //! Start executing events
    static void start();  
    //! Stop execution of events
    static void stop();
    //! Remove all events
    static void clear();
  protected:
    //static void cancel_all(Base_Signal *s);
  private:
    typedef std::deque<Base_Event*, std::allocator< Base_Event* > >::iterator Base_Event_Iterator;
    static void _run();
    static bool keep_running;
    static Ttype t; // Current time.
    static std::priority_queue<Base_Event*, 
      std::deque<Base_Event*, std::allocator<Base_Event*> >, 
      Compare_Base_Event_Times> event_queue; // Queue for the Events.
  };

  /*!

  */
  template <class ObjectType>
    class Event : public Base_Event {
  public:
    //!
    Event(ObjectType *object_pointer, void (ObjectType::*object_function_pointer)(), const Ttype delta_time) : Base_Event(delta_time) {
      po = object_pointer;
      pm = object_function_pointer;
    }
      
      //!
      virtual ~Event(){}

      //! Execute (call) the assigned function
      virtual void exec(void) { (*po.*pm)(); }  

  private:
      void (ObjectType::*pm)(); // Pointer to class member function to be executed on event expire.
      ObjectType *po; // Pointer to object who's member function is to be executed on event expire.
  };
  
  /*!

  */
  template <class ObjectType, class DataType> class Data_Event : public Base_Event {
  public:
    //!
    Data_Event(ObjectType *object_pointer, 
	       void (ObjectType::*object_function_pointer)(DataType data), 
	       DataType data, const Ttype delta_time) : Base_Event(delta_time) {
      po = object_pointer;
      pm = object_function_pointer;
      u = data;
    }

      //!
      virtual ~Data_Event(){}

      //! Execute (call) the assigned function with user data.
      virtual void exec(void) {
	(*po.*pm)(u);
      }

  private:
      void (ObjectType::*pm)(DataType data); // Pointer to class member function to be executed on event expire.
      ObjectType* po; // Pointer to object who's member function is to be executed on event expire.
      DataType u; // User data.
  };

} // namespace itpp

#endif //__events_h

