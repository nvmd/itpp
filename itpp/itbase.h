/*!
 * \file 
 * \brief Include file for the IT++ base module
 * \author Tony Ottosson
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

#ifndef ITBASE_H
#define ITBASE_H

//! \defgroup base Base Module
//@{
//! \defgroup parser Argument Parser
//! \defgroup besselfunctions Bessel Functions
//! \defgroup convertfunc Conversion Functions
//! \defgroup determinant Determinant
//! \defgroup detsource Deterministic Sources
//! \defgroup diag Diagonal Matrices and Functions
//! \defgroup errorhandlingfunc Error and Warning Handling
//! \defgroup matrix_functions Functions on Matrices
//! \defgroup hypfunc Hyperbolic Functions
//! \defgroup inverse Inverse Matrix
//! \defgroup itfile IT++ File Format
//! \defgroup logexpfunc Logarithmic and Exponential Functions
//! \defgroup matrixdecomp Matrix Decompositions
//! \defgroup miscfunc Miscellaneous Functions
//! \defgroup integration Numerical Integration
//! \defgroup optimization Numerical Optimization Routines
//! \defgroup randgen Random Number Generation
//! \defgroup reshaping Reshaping of Vectors and Matrices
//! \defgroup linearequations Solving Linear Equation Systems
//! \defgroup specmat Special Matrices
//! \defgroup statistics Statistics
//! \defgroup timers Timers
//! \defgroup trifunc Trigonometric Functions
//! \defgroup upsample Upsampling of Vectors and Matrices
//@}

#include <itpp/base/array.h>
#include <itpp/base/bessel.h>
#include <itpp/base/binary.h>
#include <itpp/base/binfile.h>
#include <itpp/base/cholesky.h>
#include <itpp/base/circular_buffer.h>
#include <itpp/base/converters.h>
#include <itpp/base/det.h>
#include <itpp/base/eigen.h>
#include <itpp/base/elmatfunc.h>
#include <itpp/base/errorfunc.h>
#include <itpp/base/factory.h>
#include <itpp/base/fastmath.h>
#include <itpp/base/gf2mat.h>
#include <itpp/base/integration.h>
#include <itpp/base/inv.h>
#include <itpp/base/itassert.h>
#include <itpp/base/itfile.h>
#include <itpp/base/itmisc.h>
#include <itpp/base/logexpfunc.h>
#include <itpp/base/ls_solve.h>
#include <itpp/base/lu.h>
#include <itpp/base/mat.h>
#include <itpp/base/matfunc.h>
#include <itpp/base/newton_search.h>
#include <itpp/base/operators.h>
#include <itpp/base/parser.h>
#include <itpp/base/qr.h>
#include <itpp/base/random.h>
#include <itpp/base/schur.h>
#include <itpp/base/smat.h>
#include <itpp/base/sort.h>
#include <itpp/base/source.h>
#include <itpp/base/specmat.h>
#include <itpp/base/stack.h>
#include <itpp/base/stat.h>
#include <itpp/base/svd.h>
#include <itpp/base/timing.h>
#include <itpp/base/trihypfunc.h>
#include <itpp/base/vec.h>

#endif // #ifndef ITBASE_H
