# File:   FindFFT.cmake
# Brief:  Find the FFT includes and library
# Author: Bogdan Cristea and Jed Brown
#
# Usage: Use the following line in your CMakeLists file: find_package(FFT)
#
# This script tries to find the path to FFT libraries and header files. The
# following variables are set:
#  FFT_VENDOR       - Possible values: ACML, MKL, FFTW3
#  FFT_INCLUDES     - where to find header file
#  FFT_LIBRARIES    - List of libraries 
#  FFT_FOUND        - True if FFT found.
#
# The original version of this script can be found at:
# https://github.com/jedbrown/cmake-modules/blob/master/FindFFTW.cmake
#
# -------------------------------------------------------------------------
#
# Copyright (C) 1995-2012  (see AUTHORS file for a list of contributors)
#
# This file is part of IT++ - a C++ library of mathematical, signal
# processing, speech processing, and communications classes and functions.
#
# IT++ is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# IT++ is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along
# with IT++.  If not, see <http://www.gnu.org/licenses/>.
#
# -------------------------------------------------------------------------

if (FFT_INCLUDES)
  # Already in cache, be silent
  set (FFT_FIND_QUIETLY TRUE)
else(FFT_INCLUDES)
  set (FFT_FOUND FALSE)
endif (FFT_INCLUDES)

#check for FFTW3
if (NOT FFT_FOUND)
  if (BLA_VENDOR MATCHES "^ACML")
    find_package(BLAS)
    if (BLAS_LIBRARIES)
      set (FFT_VENDOR "ACML")
      set (FFT_LIBRARIES ${BLAS_LIBRARIES})
      get_filename_component (PARENT_PATH ${BLAS_LIBRARIES} PATH )
      find_path (FFT_INCLUDES acml.h PATH "${PARENT_PATH}/../include")
      set (FFT_INCLUDES ${BLAS_INCLUDES})
    endif()
  elseif (BLA_VENDOR MATCHES "^Intel")
    find_package(BLAS)
    if (BLAS_LIBRARIES)
      set (FFT_VENDOR "Intel")
      set (FFT_LIBRARIES ${BLAS_LIBRARIES})
      foreach (lib ${BLAS_LIBRARIES})
        get_filename_component (PARENT_PATH ${lib} PATH )
        find_path (FFT_INCLUDES mkl_dfti.h PATH "${PARENT_PATH}/../../include")
        if (FFT_INCLUDES)
          break()
        endif()
      endforeach()
      if (FFT_INCLUDES)
        find_path(TEMP_INC fftw3.h PATH ${FFT_INCLUDES}/fftw NO_DEFAULT_PATH)
        if (TEMP_INC)
          set(FFT_INCLUDES ${TEMP_INC} ${FFT_INCLUDES})
        else()
          set(FFT_INCLUDES "")
        endif()
      endif()
      set(DFT_SEARCH_LIBS "mkl_cdft_core_dll; mkl_cdft_core")
      foreach(search_library ${DFT_SEARCH_LIBS})
        foreach (lib ${BLAS_LIBRARIES})
          get_filename_component (PARENT_PATH ${lib} PATH )
          find_library(DFT_LIBRARY
            NAMES ${search_library}
            PATHS ${PARENT_PATH})
          mark_as_advanced(DFT_LIBRARY)
          if (DFT_LIBRARY)
            set(FFT_LIBRARIES ${FFT_LIBRARIES} ${DFT_LIBRARY})
            break()
          endif()
        endforeach()
      endforeach()
      if (NOT DFT_LIBRARY)
        set(FFT_LIBRARIES "")
      endif()
    endif()
  else()
    #check for generic FFT
    find_path (FFT_INCLUDES fftw3.h)
    find_library (FFT_LIBRARIES NAMES fftw3)
    if (FFT_INCLUDES AND FFT_LIBRARIES)
      set (FFT_VENDOR "FFTW3")
    endif()
  endif()
endif()

# handle the QUIETLY and REQUIRED arguments and set FFT_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (FFT DEFAULT_MSG FFT_LIBRARIES FFT_INCLUDES)

if (FFT_LIBRARIES MATCHES "NOTFOUND$")
  set(FFT_LIBRARIES "")
endif()
if (FFT_INCLUDES MATCHES "NOTFOUND$")
  set(FFT_INCLUDES "")
endif()

mark_as_advanced (FFT_LIBRARIES FFT_INCLUDES)
