/*!
 * \file
 * \brief Templated Vector Class Implementation
 * \author Tony Ottosson, Tobias Ringstrom and Adam Piatyszek
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

#include <itpp/base/vec.h>

#if defined (HAVE_CBLAS)
#  include <itpp/base/cblas.h>
#endif


namespace itpp {

  template<>
  bool Vec<double>::set(const char *values)
  {
    std::istringstream buffer(values);
    double b, c;
    bool b_parsed = false;
    bool c_parsed = false;
    bool comma = true;
    int pos = 0, maxpos = 10;
    
    alloc(maxpos);
    
    while (buffer.peek()!=EOF) {

      switch (buffer.peek()) {
      case ',':
	it_assert(!comma, "Vec<double>::set(): Improper data string");
	buffer.seekg(1, std::ios_base::cur);
	comma = true;
	break;

      case ' ': case '\t':
	buffer.seekg(1, std::ios_base::cur);
	break;

      case ':': // reads format a:b:c or a:b
	it_assert(pos == 1, "Vec<double>::set(): Improper data string (a:b)");
	buffer.seekg(1, std::ios_base::cur);
	// parse b
	while (buffer.peek() != EOF) {
	  switch (buffer.peek()) {
	  case ' ': case '\t':
	    buffer.seekg(1, std::ios_base::cur);
	    break;

	  case ':':
	    it_assert(b_parsed, "Vec<double>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    // parse c
	    while (buffer.peek() != EOF) {
	      switch (buffer.peek()) {
	      case ' ': case '\t':
		buffer.seekg(1, std::ios_base::cur);
		break;

	      default:
		it_assert(!c_parsed, "Vec<double>::set(): Improper data string (a:b:c)");
		buffer.clear();
		buffer >> c;
		it_assert(!buffer.fail(), "Vec<double>::set(): Stream operation failed (buffer >> c)");
		c_parsed = true;
		comma = false;
	      }
	    }
	    it_assert(c_parsed, "Vec<double>::set(): Improper data string (a:b:c)");
	    break;
	    
	  default:
	    it_assert(!b_parsed, "Vec<double>::set(): Improper data string (a:b)");
	    buffer.clear();
	    buffer >> b;
	    it_assert(!buffer.fail(), "Vec<double>::set(): Stream operation failed (buffer >> b)");
	    b_parsed = true;
	    comma = false;
	  }
	}
	it_assert(b_parsed, "Vec<double>::set(): Improper data string (a:b)");

	if (c_parsed) {
	  if (b > 0 && c >= data[pos-1]) {
	    while (data[pos-1] + b <= c) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + b;
	    }
	  }
	  else if (b < 0 && c <= data[pos-1]) {
	    while (data[pos-1] + b >= c) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + b;
	    }
	  }
	  else if (b == 0 && c == data[pos-1]) {
	    break;
	  }
	  else { 
	    it_error("Vec<double>::set(): Improper data string (a:b:c)");
	  }
	} // end if (c_parsed) 
	else if (b_parsed) {
	  if (b < data[pos-1]) {
	    while (data[pos-1] > b) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] - 1.0;
	    }
	  }
	  else {
	    while (data[pos-1] < b) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + 1.0;
	    }
	  }
	} // else if (b_parsed)
	else {
	  it_error("Vec<double>::set(): Improper data string (a:b)");
	}
	break;

      default:
	if (++pos > maxpos) {
	  maxpos <<= 1;
	  set_size(maxpos, true);
	}
	buffer >> data[pos-1];
	it_assert(!buffer.fail(), "Vec<double>::set(): Stream operation failed (buffer >> data)");
	comma = false;
	break;
      }
    }
    it_assert(!comma || (pos == 0), "Vec<double>::set(): Improper data string");

    set_size(pos, true);

    return true;
  }


  template<>
  bool Vec<std::complex<double> >::set(const char *values)
  {
    std::istringstream buffer(values);
    int pos = 0, maxpos = 10;
    bool comma = true;

    alloc(maxpos);
    
    while (buffer.peek() != EOF) {
      switch (buffer.peek()) {
      case ':':
        it_error("Vec<complex>::set(): a:b:c and a:b expressions not valid for cvec");
	break;
      case ',':
	it_assert(!comma, "Vec<complex>::set(): Improper data string");
	buffer.seekg(1, std::ios_base::cur);
	comma = true;
	break;
      case ' ': case '\t':
	buffer.seekg(1, std::ios_base::cur);
	break;
      default:
	if (++pos > maxpos) {
	  maxpos <<= 1;
	  set_size(maxpos, true);
	}
	buffer >> data[pos-1];
	it_assert(!buffer.fail(), "Vec<complex>::set(): Stream operation failed (buffer >> data)");
	comma = false;
      }
    }
    it_assert(!comma || (pos == 0), "Vec<complex>::set(): Improper data string");
    set_size(pos, true);
		
    return true;
  }


  template<>
  bool Vec<bin>::set(const char *values)
  {
    std::istringstream buffer(values);
    int pos = 0, maxpos = 10;
    bool comma = true;

    alloc(maxpos);

    while (buffer.peek() != EOF) {
      switch (buffer.peek()) {
      case ':':
        it_error("Vec<bin>::set(): a:b:c and a:b expressions not valid for bvec");
	break;
      case ',':
	it_assert(!comma, "Vec<bin>::set(): Improper data string");
	buffer.seekg(1, std::ios_base::cur);
	comma = true;
	break;
      case ' ': case '\t':
	buffer.seekg(1, std::ios_base::cur);
	break;
      default:
	if (++pos > maxpos) {
	  maxpos <<= 1;
	  set_size(maxpos, true);
	}
	buffer >> data[pos-1];
	it_assert(!buffer.fail(), "Vec<bin>::set(): Stream operation failed (buffer >> data)");
	comma = false;
      }
    }
    it_assert(!comma || (pos == 0), "Vec<bin>::set(): Improper data string");

    set_size(pos, true);

    return true;
  }


  template<>
  bool Vec<int>::set(const char *values)
  {
    std::istringstream buffer(values);
    int b, c;
    bool b_parsed = false;
    bool c_parsed = false;
    bool comma = true;
    int pos = 0;
    int maxpos = 10;
    std::streamoff offset = 0;
 
    alloc(maxpos);

    while (buffer.peek() != EOF) {
      
      switch (buffer.peek()) {
      case ' ': case '\t':
	buffer.seekg(1, std::ios_base::cur);
	break;

      case ',':
	it_assert(!comma, "Vec<int>::set(): Improper data string");
	buffer.seekg(1, std::ios_base::cur);
	comma = true;
	break;

      case '+': case '-':
	buffer.seekg(1, std::ios_base::cur);
	offset--;
	break;

	// hexadecimal number or octal number or zero
      case '0':
	buffer.seekg(1, std::ios_base::cur);
	offset--;
	if (buffer.peek() != EOF) {
	  switch (buffer.peek()) {
	    // hexadecimal number
	  case 'x': case 'X':
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    do {
	      switch (buffer.peek()) {
	      case '0': case '1': case '2': case '3': case '4': case '5':
	      case '6': case '7': case '8': case '9': case 'a': case 'A':
	      case 'b': case 'B': case 'c': case 'C': case 'd': case 'D':
	      case 'e': case 'E': case 'f': case 'F':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		break;
	      default:
		it_error("Vec<int>::set(): Improper data string (hex)");
	      }
	    } while (buffer.peek() != EOF && buffer.peek() != ',' &&
		     buffer.peek() != ' ' && buffer.peek() != '\t' && 
		     buffer.peek() != ':');
	    
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    if (++pos > maxpos) {
	      maxpos <<= 1;
	      set_size(maxpos, true);
	    }
	    buffer >> std::hex >> data[pos-1];
	    comma = false;
	    break;

	    // octal number
	  case '1': case '2': case '3': case '4':
	  case '5': case '6': case '7':
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    while (buffer.peek() != EOF && buffer.peek() != ',' &&
		   buffer.peek() != ' ' && buffer.peek() != '\t' &&
		   buffer.peek() != ':') {
	      switch (buffer.peek()) {
	      case '0': case '1': case '2': case '3':
	      case '4': case '5': case '6': case '7':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		break;
	      default:
		it_error("Vec<int>::set(): Improper data string (oct)");
	      }
	    }
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    if (++pos > maxpos) {
	      maxpos <<= 1;
	      set_size(maxpos, true);
	    }
	    buffer >> std::oct >> data[pos-1];
	    comma = false;
	    break;

	    // zero
	  case ' ': case '\t': case ',': case ':':
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    if (++pos > maxpos) {
	      maxpos <<= 1;
	      set_size(maxpos, true);
	    }
	    buffer >> std::dec >> data[pos-1];
	    comma = false;
	    break;

	  default:
	    it_error("Vec<int>::set(): Improper data string");
	  }
	} // if (buffer.peek() != EOF)
	else {
	  buffer.clear();
	  buffer.seekg(offset, std::ios_base::cur);
	  offset = 0;
	  if (++pos > maxpos) {
	    maxpos <<= 1;
	    set_size(maxpos, true);
	  }
	  buffer >> std::dec >> data[pos-1];
	  comma = false;
	}
	break;

	// decimal number
      case '1': case '2': case '3': case '4': case '5':
      case '6': case '7': case '8': case '9':
	buffer.seekg(1, std::ios_base::cur);
	offset--;
	while (buffer.peek() != EOF && buffer.peek() != ',' &&
	       buffer.peek() != ' ' && buffer.peek() != '\t' &&
	       buffer.peek() != ':') {
	  switch (buffer.peek()) {
	  case '0': case '1': case '2': case '3': case '4':
	  case '5': case '6': case '7': case '8': case '9':
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    break;
	  default:
	    it_error("Vec<int>::set(): Improper data string (dec)");
	  }
	}
	buffer.clear();
	buffer.seekg(offset, std::ios_base::cur);
	offset = 0;
	if (++pos > maxpos) {
	  maxpos <<= 1;
	  set_size(maxpos, true);
	}
	buffer >> std::dec >> data[pos-1];
	comma = false;
	break;
      
	// parse format a:b:c or a:b
      case ':':
	it_assert(pos == 1, "Vec<int>::set(): Improper data string (a:b)");
	buffer.seekg(1, std::ios_base::cur);
	// parse b
	while (buffer.peek() != EOF) {
	  switch (buffer.peek()) {
	  case ' ': case '\t':
	    buffer.seekg(1, std::ios_base::cur);
	    break;
	  
	  case '+': case '-':
	    it_assert(!b_parsed, "Vec<int>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    break;
	  
	    // hexadecimal number or octal number or zero
	  case '0':
	    it_assert(!b_parsed, "Vec<int>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    if (buffer.peek() != EOF) {
	      switch (buffer.peek()) {
		// hexadecimal number
	      case 'x': case 'X':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		do {
		  switch (buffer.peek()) {
		  case '0': case '1': case '2': case '3': case '4': case '5':
		  case '6': case '7': case '8': case '9': case 'a': case 'A':
		  case 'b': case 'B': case 'c': case 'C': case 'd': case 'D':
		  case 'e': case 'E': case 'f': case 'F':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    break;
		  default:
		    it_error("Vec<int>::set(): Improper data string (a:b, hex)");
		  }
		} while (buffer.peek() != EOF && buffer.peek() != ' ' &&
			 buffer.peek() != '\t' && buffer.peek() != ':');
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::hex >> b;
		b_parsed = true;
		comma = false;
		break;

		// octal number
	      case '1': case '2': case '3': case '4':
	      case '5': case '6': case '7':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		while (buffer.peek() != EOF && buffer.peek() != ' '
		       && buffer.peek() != '\t' && buffer.peek() != ':') {
		  switch (buffer.peek()) {
		  case '0': case '1': case '2': case '3':
		  case '4': case '5': case '6': case '7':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    break;
		  default:
		    it_error("Vec<int>::set(): Improper data string (a:b, oct)");
		  }
		}
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::oct >> b;
		b_parsed = true;
		comma = false;
		break;

		// zero
	      case ' ': case '\t': case ':':
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::dec >> b;
		b_parsed = true;
		comma = false;
		break;

	      default:
		it_error("Vec<int>::set(): Improper data string (a:b)");
	      } // switch (buffer.peek())
	    } // if (buffer.peek() != EOF)
	    else {
	      buffer.clear();
	      buffer.seekg(offset, std::ios_base::cur);
	      offset = 0;
	      buffer >> std::dec >> b;
	      b_parsed = true;
	      comma = false;
	    }
	    break;
	  
	    // decimal number
	  case '1': case '2': case '3': case '4': case '5':
	  case '6': case '7': case '8': case '9':
	    it_assert(!b_parsed, "Vec<int>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    while (buffer.peek() != EOF && buffer.peek() != ' ' &&
		   buffer.peek() != '\t' && buffer.peek() != ':') {
	      switch (buffer.peek()) {
	      case '0': case '1': case '2': case '3': case '4':
	      case '5': case '6': case '7': case '8': case '9':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		break;
	      default:
		it_error("Vec<int>::set(): Improper data string (a:b, dec)");
	      }
	    }
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    buffer >> std::dec >> b;
	    b_parsed = true;
	    comma = false;
	    break;

	  case ':':
	    it_assert(b_parsed, "Vec<int>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    // parse c
	    while (buffer.peek() != EOF) {
	      switch (buffer.peek()) {
	      case ' ': case '\t':
		buffer.seekg(1, std::ios_base::cur);
		break;
	  
	      case '+': case '-':
		it_assert(!c_parsed, "Vec<int>::set(): Improper data string (a:b:c)");
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		break;
	  
		// hexadecimal number or octal number or zero
	      case '0':
		it_assert(!c_parsed, "Vec<int>::set(): Improper data string (a:b:c)");
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		if (buffer.peek() != EOF) {
		  switch (buffer.peek()) {
		    // hexadecimal number
		  case 'x': case 'X':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    do {
		      switch (buffer.peek()) {
		      case '0': case '1': case '2': case '3': case '4':
		      case '5': case '6': case '7': case '8': case '9':
		      case 'a': case 'A': case 'b': case 'B': case 'c':
		      case 'C': case 'd': case 'D': case 'e': case 'E':
		      case 'f': case 'F':
			buffer.seekg(1, std::ios_base::cur);
			offset--;
			break;
		      default:
			it_error("Vec<int>::set(): Improper data string (a:b:c, hex)");
		      }
		    } while (buffer.peek() != EOF && buffer.peek() != ' ' &&
			     buffer.peek() != '\t');
		    buffer.clear();
		    buffer.seekg(offset, std::ios_base::cur);
		    offset = 0;
		    buffer >> std::hex >> c;
		    c_parsed = true;
		    comma = false;
		    break;

		    // octal number
		  case '1': case '2': case '3': case '4':
		  case '5': case '6': case '7':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    while (buffer.peek() != EOF && buffer.peek() != ' '
			   && buffer.peek() != '\t') {
		      switch (buffer.peek()) {
		      case '0': case '1': case '2': case '3':
		      case '4': case '5': case '6': case '7':
			buffer.seekg(1, std::ios_base::cur);
			offset--;
			break;
		      default:
			it_error("Vec<int>::set(): Improper data string (a:b:c, oct)");
		      }
		    }
		    buffer.clear();
		    buffer.seekg(offset, std::ios_base::cur);
		    offset = 0;
		    buffer >> std::oct >> c;
		    c_parsed = true;
		    comma = false;
		    break;

		    // zero
		  case ' ': case '\t':
		    buffer.clear();
		    buffer.seekg(offset, std::ios_base::cur);
		    offset = 0;
		    buffer >> std::dec >> c;
		    c_parsed = true;
		    comma = false;
		    break;

		  default:
		    it_error("Vec<int>::set(): Improper data string (a:b:c)");
		  }
		}
		else {
		  buffer.clear();
		  buffer.seekg(offset, std::ios_base::cur);
		  offset = 0;
		  buffer >> std::dec >> c;
		  c_parsed = true;
		  comma = false;
		}
		break;
		
		// decimal number
	      case '1': case '2': case '3': case '4': case '5':
	      case '6': case '7': case '8': case '9':
		it_assert(!c_parsed, "Vec<int>::set(): Improper data string (a:b:c)");
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		while (buffer.peek() != EOF && buffer.peek() != ' ' &&
		       buffer.peek() != '\t') {
		  switch (buffer.peek()) {
		  case '0': case '1': case '2': case '3': case '4':
		  case '5': case '6': case '7': case '8': case '9':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    break;
		  default:
		    it_error("Vec<int>::set(): Improper data string (a:b:c, dec)");
		  }
		}
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::dec >> c;
		c_parsed = true;
		comma = false;
		break;

	      default:
		it_error("Vec<int>::set(): Improper data string (a:b:c)");
	      } // switch (buffer.peek())
	    } // while (buffer.peek() != EOF)
	    it_assert(c_parsed, "Vec<int>::set(): Improper data string (a:b:c)");
	    break;
      
	  default:
	    it_error("Vec<int>::set(): Improper data string (a:b)");
	  } // switch (buffer.peek())
	} // while (buffer.peek() != EOF)

	if (c_parsed) {
	  if (b > 0 && c >= data[pos-1]) {
	    while (data[pos-1] + b <= c) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + b;
	    }
	  }
	  else if (b < 0 && c <= data[pos-1]) {
	    while (data[pos-1] + b >= c) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + b;
	    }
	  }
	  else if (b == 0 && c == data[pos-1]) {
	    break;
	  }
	  else { 
	    it_error("Vec<int>::set(): Improper data string (a:b:c)");
	  }
	} // if (c_parsed) 
	else if (b_parsed) {
	  if (b < data[pos-1]) {
	    while (data[pos-1] > b) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] - 1;
	    }
	  }
	  else {
	    while (data[pos-1] < b) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + 1;
	    }
	  }
	} // else if (b_parsed)
	else {
	  it_error("Vec<int>::set(): Improper data string (a:b)");
	}
	break;
      
      default:
	it_error("Vec<int>::set(): Improper data string");
      }
    }
    it_assert(!comma || (pos == 0), "Vec<int>::set(): Improper data string");

    // resize the parsed vector to its final length
    set_size(pos, true);

    return true;
  }

  template<>
  bool Vec<short int>::set(const char *values)
  {
    std::istringstream buffer(values);
    short int b, c;
    bool b_parsed = false;
    bool c_parsed = false;
    bool comma = true;
    int pos = 0;
    int maxpos = 10;
    std::streamoff offset = 0;
 
    alloc(maxpos);

    while (buffer.peek() != EOF) {
      
      switch (buffer.peek()) {
      case ' ': case '\t':
	buffer.seekg(1, std::ios_base::cur);
	break;

      case ',':
	it_assert(!comma, "Vec<short>::set(): Improper data string");
	buffer.seekg(1, std::ios_base::cur);
	comma = true;
	break;

      case '+': case '-':
	buffer.seekg(1, std::ios_base::cur);
	offset--;
	break;

	// hexadecimal number or octal number or zero
      case '0':
	buffer.seekg(1, std::ios_base::cur);
	offset--;
	if (buffer.peek() != EOF) {
	  switch (buffer.peek()) {
	    // hexadecimal number
	  case 'x': case 'X':
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    do {
	      switch (buffer.peek()) {
	      case '0': case '1': case '2': case '3': case '4': case '5':
	      case '6': case '7': case '8': case '9': case 'a': case 'A':
	      case 'b': case 'B': case 'c': case 'C': case 'd': case 'D':
	      case 'e': case 'E': case 'f': case 'F':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		break;
	      default:
		it_error("Vec<short>::set(): Improper data string (hex)");
	      }
	    } while (buffer.peek() != EOF && buffer.peek() != ',' &&
		     buffer.peek() != ' ' && buffer.peek() != '\t' && 
		     buffer.peek() != ':');
	    
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    if (++pos > maxpos) {
	      maxpos <<= 1;
	      set_size(maxpos, true);
	    }
	    buffer >> std::hex >> data[pos-1];
	    comma = false;
	    break;

	    // octal number
	  case '1': case '2': case '3': case '4':
	  case '5': case '6': case '7':
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    while (buffer.peek() != EOF && buffer.peek() != ',' &&
		   buffer.peek() != ' ' && buffer.peek() != '\t' &&
		   buffer.peek() != ':') {
	      switch (buffer.peek()) {
	      case '0': case '1': case '2': case '3':
	      case '4': case '5': case '6': case '7':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		break;
	      default:
		it_error("Vec<short>::set(): Improper data string (oct)");
	      }
	    }
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    if (++pos > maxpos) {
	      maxpos <<= 1;
	      set_size(maxpos, true);
	    }
	    buffer >> std::oct >> data[pos-1];
	    comma = false;
	    break;

	    // zero
	  case ' ': case '\t': case ',': case ':':
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    if (++pos > maxpos) {
	      maxpos <<= 1;
	      set_size(maxpos, true);
	    }
	    buffer >> std::dec >> data[pos-1];
	    comma = false;
	    break;

	  default:
	    it_error("Vec<short>::set(): Improper data string");
	  }
	} // if (buffer.peek() != EOF)
	else {
	  buffer.clear();
	  buffer.seekg(offset, std::ios_base::cur);
	  offset = 0;
	  if (++pos > maxpos) {
	    maxpos <<= 1;
	    set_size(maxpos, true);
	  }
	  buffer >> std::dec >> data[pos-1];
	  comma = false;
	}
	break;

	// decimal number
      case '1': case '2': case '3': case '4': case '5':
      case '6': case '7': case '8': case '9':
	buffer.seekg(1, std::ios_base::cur);
	offset--;
	while (buffer.peek() != EOF && buffer.peek() != ',' &&
	       buffer.peek() != ' ' && buffer.peek() != '\t' &&
	       buffer.peek() != ':') {
	  switch (buffer.peek()) {
	  case '0': case '1': case '2': case '3': case '4':
	  case '5': case '6': case '7': case '8': case '9':
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    break;
	  default:
	    it_error("Vec<short>::set(): Improper data string (dec)");
	  }
	}
	buffer.clear();
	buffer.seekg(offset, std::ios_base::cur);
	offset = 0;
	if (++pos > maxpos) {
	  maxpos <<= 1;
	  set_size(maxpos, true);
	}
	buffer >> std::dec >> data[pos-1];
	comma = false;
	break;
      
	// parse format a:b:c or a:b
      case ':':
	it_assert(pos == 1, "Vec<short>::set(): Improper data string (a:b)");
	buffer.seekg(1, std::ios_base::cur);
	// parse b
	while (buffer.peek() != EOF) {
	  switch (buffer.peek()) {
	  case ' ': case '\t':
	    buffer.seekg(1, std::ios_base::cur);
	    break;
	  
	  case '+': case '-':
	    it_assert(!b_parsed, "Vec<short>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    break;
	  
	    // hexadecimal number or octal number or zero
	  case '0':
	    it_assert(!b_parsed, "Vec<short>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    if (buffer.peek() != EOF) {
	      switch (buffer.peek()) {
		// hexadecimal number
	      case 'x': case 'X':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		do {
		  switch (buffer.peek()) {
		  case '0': case '1': case '2': case '3': case '4': case '5':
		  case '6': case '7': case '8': case '9': case 'a': case 'A':
		  case 'b': case 'B': case 'c': case 'C': case 'd': case 'D':
		  case 'e': case 'E': case 'f': case 'F':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    break;
		  default:
		    it_error("Vec<short>::set(): Improper data string (a:b, hex)");
		  }
		} while (buffer.peek() != EOF && buffer.peek() != ' ' &&
			 buffer.peek() != '\t' && buffer.peek() != ':');
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::hex >> b;
		b_parsed = true;
		comma = false;
		break;

		// octal number
	      case '1': case '2': case '3': case '4':
	      case '5': case '6': case '7':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		while (buffer.peek() != EOF && buffer.peek() != ' '
		       && buffer.peek() != '\t' && buffer.peek() != ':') {
		  switch (buffer.peek()) {
		  case '0': case '1': case '2': case '3':
		  case '4': case '5': case '6': case '7':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    break;
		  default:
		    it_error("Vec<short>::set(): Improper data string (a:b, oct)");
		  }
		}
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::oct >> b;
		b_parsed = true;
		comma = false;
		break;

		// zero
	      case ' ': case '\t': case ':':
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::dec >> b;
		b_parsed = true;
		comma = false;
		break;

	      default:
		it_error("Vec<short>::set(): Improper data string (a:b)");
	      } // switch (buffer.peek())
	    } // if (buffer.peek() != EOF)
	    else {
	      buffer.clear();
	      buffer.seekg(offset, std::ios_base::cur);
	      offset = 0;
	      buffer >> std::dec >> b;
	      b_parsed = true;
	      comma = false;
	    }
	    break;
	  
	    // decimal number
	  case '1': case '2': case '3': case '4': case '5':
	  case '6': case '7': case '8': case '9':
	    it_assert(!b_parsed, "Vec<short>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    offset--;
	    while (buffer.peek() != EOF && buffer.peek() != ' ' &&
		   buffer.peek() != '\t' && buffer.peek() != ':') {
	      switch (buffer.peek()) {
	      case '0': case '1': case '2': case '3': case '4':
	      case '5': case '6': case '7': case '8': case '9':
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		break;
	      default:
		it_error("Vec<short>::set(): Improper data string (a:b, dec)");
	      }
	    }
	    buffer.clear();
	    buffer.seekg(offset, std::ios_base::cur);
	    offset = 0;
	    buffer >> std::dec >> b;
	    b_parsed = true;
	    comma = false;
	    break;

	  case ':':
	    it_assert(b_parsed, "Vec<short>::set(): Improper data string (a:b)");
	    buffer.seekg(1, std::ios_base::cur);
	    // parse c
	    while (buffer.peek() != EOF) {
	      switch (buffer.peek()) {
	      case ' ': case '\t':
		buffer.seekg(1, std::ios_base::cur);
		break;
	  
	      case '+': case '-':
		it_assert(!c_parsed, "Vec<short>::set(): Improper data string (a:b:c)");
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		break;
	  
		// hexadecimal number or octal number or zero
	      case '0':
		it_assert(!c_parsed, "Vec<short>::set(): Improper data string (a:b:c)");
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		if (buffer.peek() != EOF) {
		  switch (buffer.peek()) {
		    // hexadecimal number
		  case 'x': case 'X':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    do {
		      switch (buffer.peek()) {
		      case '0': case '1': case '2': case '3': case '4':
		      case '5': case '6': case '7': case '8': case '9':
		      case 'a': case 'A': case 'b': case 'B': case 'c':
		      case 'C': case 'd': case 'D': case 'e': case 'E':
		      case 'f': case 'F':
			buffer.seekg(1, std::ios_base::cur);
			offset--;
			break;
		      default:
			it_error("Vec<short>::set(): Improper data string (a:b:c, hex)");
		      }
		    } while (buffer.peek() != EOF && buffer.peek() != ' ' &&
			     buffer.peek() != '\t');
		    buffer.clear();
		    buffer.seekg(offset, std::ios_base::cur);
		    offset = 0;
		    buffer >> std::hex >> c;
		    c_parsed = true;
		    comma = false;
		    break;

		    // octal number
		  case '1': case '2': case '3': case '4':
		  case '5': case '6': case '7':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    while (buffer.peek() != EOF && buffer.peek() != ' '
			   && buffer.peek() != '\t') {
		      switch (buffer.peek()) {
		      case '0': case '1': case '2': case '3':
		      case '4': case '5': case '6': case '7':
			buffer.seekg(1, std::ios_base::cur);
			offset--;
			break;
		      default:
			it_error("Vec<short>::set(): Improper data string (a:b:c, oct)");
		      }
		    }
		    buffer.clear();
		    buffer.seekg(offset, std::ios_base::cur);
		    offset = 0;
		    buffer >> std::oct >> c;
		    c_parsed = true;
		    comma = false;
		    break;

		    // zero
		  case ' ': case '\t':
		    buffer.clear();
		    buffer.seekg(offset, std::ios_base::cur);
		    offset = 0;
		    buffer >> std::dec >> c;
		    c_parsed = true;
		    comma = false;
		    break;

		  default:
		    it_error("Vec<short>::set(): Improper data string (a:b:c)");
		  }
		}
		else {
		  buffer.clear();
		  buffer.seekg(offset, std::ios_base::cur);
		  offset = 0;
		  buffer >> std::dec >> c;
		  c_parsed = true;
		  comma = false;
		}
		break;
		
		// decimal number
	      case '1': case '2': case '3': case '4': case '5':
	      case '6': case '7': case '8': case '9':
		it_assert(!c_parsed, "Vec<short>::set(): Improper data string (a:b:c)");
		buffer.seekg(1, std::ios_base::cur);
		offset--;
		while (buffer.peek() != EOF && buffer.peek() != ' ' &&
		       buffer.peek() != '\t') {
		  switch (buffer.peek()) {
		  case '0': case '1': case '2': case '3': case '4':
		  case '5': case '6': case '7': case '8': case '9':
		    buffer.seekg(1, std::ios_base::cur);
		    offset--;
		    break;
		  default:
		    it_error("Vec<short>::set(): Improper data string (a:b:c, dec)");
		  }
		}
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::dec >> c;
		c_parsed = true;
		comma = false;
		break;

	      default:
		it_error("Vec<short>::set(): Improper data string (a:b:c)");
	      } // switch (buffer.peek())
	    } // while (buffer.peek() != EOF)
	    it_assert(c_parsed, "Vec<short>::set(): Improper data string (a:b:c)");
	    break;
      
	  default:
	    it_error("Vec<short>::set(): Improper data string (a:b)");
	  } // switch (buffer.peek())
	} // while (buffer.peek() != EOF)

	if (c_parsed) {
	  if (b > 0 && c >= data[pos-1]) {
	    while (data[pos-1] + b <= c) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + b;
	    }
	  }
	  else if (b < 0 && c <= data[pos-1]) {
	    while (data[pos-1] + b >= c) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + b;
	    }
	  }
	  else if (b == 0 && c == data[pos-1]) {
	    break;
	  }
	  else { 
	    it_error("Vec<short>::set(): Improper data string (a:b:c)");
	  }
	} // if (c_parsed) 
	else if (b_parsed) {
	  if (b < data[pos-1]) {
	    while (data[pos-1] > b) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] - 1;
	    }
	  }
	  else {
	    while (data[pos-1] < b) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + 1;
	    }
	  }
	} // else if (b_parsed)
	else {
	  it_error("Vec<short>::set(): Improper data string (a:b)");
	}
	break;
      
      default:
	it_error("Vec<short>::set(): Improper data string");
      }
    }
    it_assert(!comma || (pos == 0), "Vec<short>::set(): Improper data string");

    // resize the parsed vector to its final length
    set_size(pos, true);

    return true;
  }


#if defined(HAVE_CBLAS)
  template<>
  double dot(const vec &v1, const vec &v2)
  {
    it_assert1(v1.datasize == v2.datasize, "vec::dot: wrong sizes");
    double r=0.0;

    r= cblas_ddot(v1.datasize, v1.data, 1, v2.data, 1);

    return r;
  }
 
  template<>
  std::complex<double> dot(const cvec &v1, const cvec &v2)
  {
    it_assert1(v1.datasize == v2.datasize, "cvec::dot: wrong sizes");
    std::complex<double> r=0.0;
		
    cblas_zdotu_sub(v1.datasize, v1.data, 1, v2.data, 1, &r);

    return r;
  }

#endif // HAVE_CBLAS

  template<> 
  bvec Vec<std::complex<double> >::operator==(const std::complex<double>) const
  { 
    it_error("operator==: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec Vec<std::complex<double> >::operator!=(const std::complex<double>) const
  { 
    it_error("operator!=: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec Vec<std::complex<double> >::operator<=(const std::complex<double>) const
  { 
    it_error("operator<=: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec Vec<std::complex<double> >::operator>(const std::complex<double>) const
  { 
    it_error("operator>: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec Vec<std::complex<double> >::operator<(const std::complex<double>) const
  { 
    it_error("operator<: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<> 
  bvec Vec<std::complex<double> >::operator>=(const std::complex<double>) const
  { 
    it_error("operator>=: not implemented for complex");
    bvec temp;
    return temp;
  }

  template<>
  Mat<std::complex<double> > Vec<std::complex<double> >::hermitian_transpose() const
  {
    Mat<std::complex<double> > temp(1, datasize);
    for (int i=0; i<datasize; i++)
      temp(i) = std::conj(data[i]);

    return temp;
  }


  //---------------------------------------------------------------------
  // Instantiations
  //---------------------------------------------------------------------

  //--------- class instantiations -------------

  template class Vec<double>;
  template class Vec<int>;
  template class Vec<short int>;
  template class Vec<std::complex<double> >;
  template class Vec<bin>;

  //------------- Addition operator ----------

  template const vec operator+(const vec &v1, const vec &v2);
  template const cvec operator+(const cvec &v1, const cvec &v2);
  template const ivec operator+(const ivec &v1, const ivec &v2);
  template const svec operator+(const svec &v1, const svec &v2);
  template const bvec operator+(const bvec &v1, const bvec &v2);

  template const vec operator+(const vec &v1, double t);
  template const cvec operator+(const cvec &v1, std::complex<double> t);
  template const ivec operator+(const ivec &v1, int t);
  template const svec operator+(const svec &v1, short t);
  template const bvec operator+(const bvec &v1, bin t);

  template const vec operator+(double t, const vec &v1);
  template const cvec operator+(std::complex<double> t, const cvec &v1);
  template const ivec operator+(int t, const ivec &v1);
  template const svec operator+(short t, const svec &v1);
  template const bvec operator+(bin t, const bvec &v1);

  //------------- Subraction operator ----------

  template const vec operator-(const vec &v1, const vec &v2);
  template const cvec operator-(const cvec &v1, const cvec &v2);
  template const ivec operator-(const ivec &v1, const ivec &v2);
  template const svec operator-(const svec &v1, const svec &v2);
  template const bvec operator-(const bvec &v1, const bvec &v2);

  template const vec operator-(const vec &v, double t);
  template const cvec operator-(const cvec &v, std::complex<double> t);
  template const ivec operator-(const ivec &v, int t);
  template const svec operator-(const svec &v, short t);
  template const bvec operator-(const bvec &v, bin t);

  template const vec operator-(double t, const vec &v);
  template const cvec operator-(std::complex<double> t, const cvec &v);
  template const ivec operator-(int t, const ivec &v);
  template const svec operator-(short t, const svec &v);
  template const bvec operator-(bin t, const bvec &v);

  //---------- Unary minus -------------

  template const vec operator-(const vec &v);
  template const cvec operator-(const cvec &v);
  template const ivec operator-(const ivec &v);
  template const svec operator-(const svec &v);
  template const bvec operator-(const bvec &v);

  //------------- Multiplication operator ----------

#if !defined(HAVE_CBLAS)
  template double dot(const vec &v1, const vec &v2);
  template std::complex<double> dot(const cvec &v1, const cvec &v2);
#endif
  template int dot(const ivec &v1, const ivec &v2);
  template short dot(const svec &v1, const svec &v2);
  template bin dot(const bvec &v1, const bvec &v2);

  template int operator*(const ivec &v1, const ivec &v2);
  template short operator*(const svec &v1, const svec &v2);
  template bin operator*(const bvec &v1, const bvec &v2);

  template const mat outer_product(const vec &v1, const vec &v2);
  template const cmat outer_product(const cvec &v1, const cvec &v2);
  template const imat outer_product(const ivec &v1, const ivec &v2);
  template const smat outer_product(const svec &v1, const svec &v2);
  template const bmat outer_product(const bvec &v1, const bvec &v2);

  template const vec operator*(const vec &v, double t);
  template const cvec operator*(const cvec &v, std::complex<double> t);
  template const ivec operator*(const ivec &v, int t);
  template const svec operator*(const svec &v, short t);
  template const bvec operator*(const bvec &v, bin t);

  template const vec operator*(double t, const vec &v);
  template const cvec operator*(std::complex<double> t, const cvec &v);
  template const ivec operator*(int t, const ivec &v);
  template const svec operator*(short t, const svec &v);
  template const bvec operator*(bin t, const bvec &v);

  //------------- Elementwise Multiplication operator (two vectors) ----------

  template const vec elem_mult(const vec &v1, const vec &v2);
  template const cvec elem_mult(const cvec &v1, const cvec &v2);
  template const ivec elem_mult(const ivec &v1, const ivec &v2);
  template const svec elem_mult(const svec &v1, const svec &v2);
  template const bvec elem_mult(const bvec &v1, const bvec &v2);

  //------------- Elementwise Multiplication operator (three vectors) ----------

  template const vec elem_mult(const vec &v1, const vec &v2, const vec &v3);
  template const cvec elem_mult(const cvec &v1, const cvec &v2, const cvec &v3);
  template const ivec elem_mult(const ivec &v1, const ivec &v2, const ivec &v3);
  template const svec elem_mult(const svec &v1, const svec &v2, const svec &v3);
  template const bvec elem_mult(const bvec &v1, const bvec &v2, const bvec &v3);

  //------------- Elementwise Multiplication operator (four vectors) ----------

  template const vec elem_mult(const vec &v1, const vec &v2, const vec &v3, const vec &v4);
  template const cvec elem_mult(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4);
  template const ivec elem_mult(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4);
  template const svec elem_mult(const svec &v1, const svec &v2, const svec &v3, const svec &v4);
  template const bvec elem_mult(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4);

  //------------- Division operator ----------

  template const vec operator/(const vec &v, double t);
  template const cvec operator/(const cvec &v, std::complex<double> t);
  template const ivec operator/(const ivec &v, int t);
  template const svec operator/(const svec &v, short t);
  template const bvec operator/(const bvec &v, bin t);

  template const vec operator/(const double t, const vec &v);
  template const cvec operator/(const std::complex<double> t, const cvec &v);
  template const ivec operator/(const int t, const ivec &v);
  template const svec operator/(const short t, const svec &v);
  template const bvec operator/(const bin t, const bvec &v);

  //------------- Elementwise Division operator ----------

  template const vec elem_div(const vec &v1, const vec &v2);
  template const cvec elem_div(const cvec &v1, const cvec &v2);
  template const ivec elem_div(const ivec &v1, const ivec &v2);
  template const svec elem_div(const svec &v1, const svec &v2);
  template const bvec elem_div(const bvec &v1, const bvec &v2);

  template const vec elem_div(const double t, const vec &v);
  template const cvec elem_div(const std::complex<double> t, const cvec &v);
  template const ivec elem_div(const int t, const ivec &v);
  template const svec elem_div(const short t, const svec &v);
  template const bvec elem_div(const bin t, const bvec &v);

  //--------------------- concat operator -----------------

  template const vec concat(const vec &v, const double a);
  template const cvec concat(const cvec &v, const std::complex<double> a);
  template const ivec concat(const ivec &v, const int a);
  template const svec concat(const svec &v, const short a);
  template const bvec concat(const bvec &v, const bin a);

  template const vec concat(const double a, const vec &v);
  template const cvec concat(const std::complex<double> a, const cvec &v);
  template const ivec concat(const int a, const ivec &v);
  template const svec concat(const short a, const svec &v);
  template const bvec concat(const bin a, const bvec &v);

  template const vec concat(const vec &v1, const vec &v2);
  template const cvec concat(const cvec &v1, const cvec &v2);
  template const ivec concat(const ivec &v1, const ivec &v2);
  template const svec concat(const svec &v1, const svec &v2);
  template const bvec concat(const bvec &v1, const bvec &v2);

  template const vec concat(const vec &v1, const vec &v2, const vec &v3);
  template const cvec concat(const cvec &v1, const cvec &v2, const cvec &v3);
  template const ivec concat(const ivec &v1, const ivec &v2, const ivec &v3);
  template const svec concat(const svec &v1, const svec &v2, const svec &v3);
  template const bvec concat(const bvec &v1, const bvec &v2, const bvec &v3);

  template const vec concat(const vec &v1, const vec &v2, const vec &v3, const vec &v4);
  template const cvec concat(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4);
  template const ivec concat(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4);
  template const svec concat(const svec &v1, const svec &v2, const svec &v3, const svec &v4);
  template const bvec concat(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4);

  template const vec concat(const vec &v1, const vec &v2, const vec &v3, const vec &v4, const vec &v5);
  template const cvec concat(const cvec &v1, const cvec &v2, const cvec &v3, const cvec &v4, const cvec &v5);
  template const ivec concat(const ivec &v1, const ivec &v2, const ivec &v3, const ivec &v4, const ivec &v5);
  template const svec concat(const svec &v1, const svec &v2, const svec &v3, const svec &v4, const svec &v5);
  template const bvec concat(const bvec &v1, const bvec &v2, const bvec &v3, const bvec &v4, const bvec &v5);

  // -------------- output stream --------------------

  template std::ostream &operator<<(std::ostream& os, const vec &vect);
  template std::ostream &operator<<(std::ostream& os, const cvec &vect);
  template std::ostream &operator<<(std::ostream& os, const svec &vect);
  template std::ostream &operator<<(std::ostream& os, const ivec &vect);
  template std::ostream &operator<<(std::ostream& os, const bvec &vect);

  // -------------- input stream --------------------

  template std::istream &operator>>(std::istream& is, vec &vect);
  template std::istream &operator>>(std::istream& is, cvec &vect);
  template std::istream &operator>>(std::istream& is, svec &vect);
  template std::istream &operator>>(std::istream& is, ivec &vect);
  template std::istream &operator>>(std::istream& is, bvec &vect);

} // namespace itpp
