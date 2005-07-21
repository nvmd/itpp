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
  \brief Include file for the it++ communication protocols library

  Developed at the Dept. of Signals and Systems, Chalmers University of Technology,
  SE-412 96 Göteborg.

  $Revision$

  $Date$
*/

#ifndef __itprotocol_h
#define __itprotocol_h

#include <itpp/itcomm.h>

#include <itpp/protocol/events.h>
#include <itpp/protocol/signals_slots.h>

#include <itpp/protocol/packet.h>
#include <itpp/protocol/packet_generator.h>
#include <itpp/protocol/packet_channel.h>

#include <itpp/protocol/front_drop_queue.h>

#include <itpp/protocol/selective_repeat.h>

#include <itpp/protocol/tcp.h>
#include <itpp/protocol/tcp_client_server.h>


#endif // __itprotocol_h
