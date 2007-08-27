/*!
 * \file
 * \brief Base class for class factories and memory allocation functions
 * \author Johan Bergman and Adam Piatyszek
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

#include <itpp/base/factory.h>

namespace itpp {

  template<>
  void create_elements(unsigned char* &ptr, int n, const Factory &)
  {
    char* data_char = new char[sizeof(unsigned char) * n];
    ptr = reinterpret_cast<unsigned char*>(data_char);
  }

  template<>
  void create_elements(bin* &ptr, int n, const Factory &)
  {
    char* data_char = new char[sizeof(bin) * n];
    ptr = reinterpret_cast<bin*>(data_char);
  }

  template<>
  void create_elements(short int* &ptr, int n, const Factory &)
  {
    char* data_char = new char[sizeof(short int) * n];
    ptr = reinterpret_cast<short int*>(data_char);
  }

  template<>
  void create_elements(int* &ptr, int n, const Factory &)
  {
    char* data_char = new char[sizeof(int) * n];
    ptr = reinterpret_cast<int*>(data_char);
  }

  // The following specialisation allocate memory using 16-byte alignment
  template<>
  void create_elements(double* &ptr, int n, const Factory &)
  {
    void *p0 = ::operator new(sizeof(double) * n + 16);
    void *p1 = reinterpret_cast<void*>((reinterpret_cast<std::size_t>(p0) + 16)
				       & (~(std::size_t(15))));
    *(reinterpret_cast<void**>(p1) - 1) = p0;
    ptr = reinterpret_cast<double*>(p1);
  }

  // The following specialisation allocate memory using 16-byte alignment
  template<>
  void create_elements(std::complex<double>* &ptr, int n, const Factory &)
  {
    void *p0 = ::operator new(sizeof(std::complex<double>) * n + 16);
    void *p1 = reinterpret_cast<void*>((reinterpret_cast<std::size_t>(p0) + 16)
				       & (~(std::size_t(15))));
    *(reinterpret_cast<void**>(p1) - 1) = p0;
    ptr = reinterpret_cast<std::complex<double>*>(p1);
  }


  template<>
  void destroy_elements(unsigned char* &ptr, int n)
  {
    if (ptr) {
      delete [] reinterpret_cast<char*>(ptr);
      ptr = 0;
    }
  }

  template<>
  void destroy_elements(bin* &ptr, int n)
  {
    if (ptr) {
      delete [] reinterpret_cast<char*>(ptr);
      ptr = 0;
    }
  }

  template<>
  void destroy_elements(short int* &ptr, int n)
  {
    if (ptr) {
      delete [] reinterpret_cast<char*>(ptr);
      ptr = 0;
    }
  }

  template<>
  void destroy_elements(int* &ptr, int n)
  {
    if (ptr) {
      delete [] reinterpret_cast<char*>(ptr);
      ptr = 0;
    }
  }


  template<>
  void destroy_elements(double* &ptr, int n)
  {
    if (ptr) {
      void *p = *(reinterpret_cast<void**>(ptr) - 1);
      ::operator delete(p);
      ptr = 0;
    }
  }

  template<>
  void destroy_elements(std::complex<double>* &ptr, int n)
  {
    if (ptr) {
      void *p = *(reinterpret_cast<void**>(ptr) - 1);
      ::operator delete(p);
      ptr = 0;
    }
  }

}
