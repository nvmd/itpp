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
  \brief Include file for the it++ communication library

  Developed at the Dept. of Signals and Systems, Chalmers University of Technology,
  SE-412 96 Göteborg.

  $Revision$

  $Date$
*/

#ifndef __itcomm_h
#define __itcomm_h

#include <itpp/itbase.h>

#include <itpp/comm/modulator.h>
#include <itpp/comm/channel.h>
#include <itpp/comm/interleave.h>
#include <itpp/comm/spread.h>
#include <itpp/comm/commfunc.h>
#include <itpp/comm/galois.h>
#include <itpp/comm/channel_code.h>
#include <itpp/comm/hammcode.h>
#include <itpp/comm/convcode.h>
#include <itpp/comm/punct_convcode.h>
#include <itpp/comm/bch.h>
#include <itpp/comm/reedsolomon.h>
#include <itpp/comm/egolay.h>
#include <itpp/comm/crc.h>
#include <itpp/comm/error_counters.h>
#include <itpp/comm/sequence.h>
#include <itpp/comm/ofdm.h>
#include <itpp/comm/rec_syst_conv_code.h>
#include <itpp/comm/turbo.h>
#include <itpp/comm/pulse_shape.h>

#endif // __itcomm_h
