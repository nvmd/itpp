/*!
 * \file 
 * \brief Include file for the IT++ signal-processing module
 * \author Adam Piatyszek
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

#ifndef ITSIGNAL_H
#define ITSIGNAL_H

//! \defgroup signal Signal Processing Module
//@{
//! \defgroup dct Discrete Cosine Transform (DCT)
//! \defgroup fft Fast Fourier Transform (FFT)
//! \defgroup fht Fast Hadamard Transform (FHT)
//! \defgroup fastica Fast Independent Component Analysis
//! \defgroup filters Filtering
//! \defgroup poly Polynomial Functions
//! \defgroup sigproc Signal Processing
//! \defgroup windfunc Windowing
//@}

#include <itpp/itbase.h>
#include <itpp/signal/fastica.h>
#include <itpp/signal/filter.h>
#include <itpp/signal/filter_design.h>
#include <itpp/signal/freq_filt.h>
#include <itpp/signal/poly.h>
#include <itpp/signal/sigfun.h>
#include <itpp/signal/transforms.h>
#include <itpp/signal/window.h>

#endif // #ifndef ITSIGNAL_H
