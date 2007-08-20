/*!
 * \file
 * \brief Templated Vector Class Implementation
 * \author Tony Ottosson, Tobias Ringstrom, Adam Piatyszek and Conrad Sanderson
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

#include <itpp/base/vec.h>


namespace itpp {


  template<class Num_T>
  std::string Vec<Num_T>::replace_commas(const std::string &str_in)
  {
    // copy an input sting into a local variable str
    std::string str(str_in);
    // find first occurence of comma in string str
    std::string::size_type index = str.find(',', 0);
    while (index != std::string::npos) {
      // replace character at position index with space
      str.replace(index, 1, 1, ' ');
      // find next occurence of comma in string str
      index = str.find(',', index);
    }
    return str;
  }


  template<>
  bool Vec<double>::set(const std::string &str)
  {
    std::istringstream buffer(replace_commas(str));
    double b, c, eps_margin;
    bool b_parsed = false;
    bool c_parsed = false;
    int pos = 0, maxpos = 10;

    free();
    alloc(maxpos);

    while (buffer.peek()!=EOF) {

      switch (buffer.peek()) {
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
	  }
	}
	it_assert(b_parsed, "Vec<double>::set(): Improper data string (a:b)");

	if (c_parsed) {
	  // Adding this margin fixes precision problems in e.g. "0:0.2:3",
	  // where the last value was 2.8 instead of 3.
	  eps_margin = std::fabs((c - data[pos-1]) / b) * eps;
	  if (b > 0 && c >= data[pos-1]) {
	    while (data[pos-1] + b <= c + eps_margin) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] + b;
	    }
	  }
	  else if (b < 0 && c <= data[pos-1]) {
	    while (data[pos-1] + b >= c - eps_margin) {
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
	  eps_margin = std::fabs(b - data[pos-1]) * eps;
	  if (b < data[pos-1]) {
	    while (data[pos-1] -1.0 >= b - eps_margin) {
	      if (++pos > maxpos) {
		maxpos <<= 1;
		set_size(maxpos, true);
	      }
	      data[pos-1] = data[pos-2] - 1.0;
	    }
	  }
	  else {
	    while (data[pos-1] + 1.0 <= b + eps_margin) {
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
	break;
      }
    }
    set_size(pos, true);

    return true;
  }


  template<>
  bool Vec<std::complex<double> >::set(const std::string &str)
  {
    std::istringstream buffer(str);
    int pos = 0, maxpos = 10;
    bool comma = true;

    free();
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
  bool Vec<bin>::set(const std::string &str)
  {
    std::istringstream buffer(replace_commas(str));
    int pos = 0, maxpos = 10;

    free();
    alloc(maxpos);

    while (buffer.peek() != EOF) {
      switch (buffer.peek()) {
      case ':':
        it_error("Vec<bin>::set(): a:b:c and a:b expressions not valid for bvec");
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
      }
    }
    set_size(pos, true);

    return true;
  }


  template<>
  bool Vec<int>::set(const std::string &str)
  {
    std::istringstream buffer(replace_commas(str));
    int b, c;
    bool b_parsed = false;
    bool c_parsed = false;
    int pos = 0;
    int maxpos = 10;
    std::streamoff offset = 0;

    free();
    alloc(maxpos);

    while (buffer.peek() != EOF) {

      switch (buffer.peek()) {
      case ' ': case '\t':
	buffer.seekg(1, std::ios_base::cur);
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
		break;

		// zero
	      case ' ': case '\t': case ':':
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::dec >> b;
		b_parsed = true;
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
		    break;

		    // zero
		  case ' ': case '\t':
		    buffer.clear();
		    buffer.seekg(offset, std::ios_base::cur);
		    offset = 0;
		    buffer >> std::dec >> c;
		    c_parsed = true;
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
    // resize the parsed vector to its final length
    set_size(pos, true);

    return true;
  }

  template<>
  bool Vec<short int>::set(const std::string &str)
  {
    std::istringstream buffer(replace_commas(str));
    short int b, c;
    bool b_parsed = false;
    bool c_parsed = false;
    int pos = 0;
    int maxpos = 10;
    std::streamoff offset = 0;

    free();
    alloc(maxpos);

    while (buffer.peek() != EOF) {

      switch (buffer.peek()) {
      case ' ': case '\t':
	buffer.seekg(1, std::ios_base::cur);
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
		break;

		// zero
	      case ' ': case '\t': case ':':
		buffer.clear();
		buffer.seekg(offset, std::ios_base::cur);
		offset = 0;
		buffer >> std::dec >> b;
		b_parsed = true;
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
		    break;

		    // zero
		  case ' ': case '\t':
		    buffer.clear();
		    buffer.seekg(offset, std::ios_base::cur);
		    offset = 0;
		    buffer >> std::dec >> c;
		    c_parsed = true;
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
    // resize the parsed vector to its final length
    set_size(pos, true);

    return true;
  }


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

#if !defined(HAVE_BLAS)
  template double dot(const vec &v1, const vec &v2);
#if defined(HAVE_MKL) && defined(_MSC_VER)
  template std::complex<double> dot(const cvec &v1, const cvec &v2);
#endif
#endif
  template int dot(const ivec &v1, const ivec &v2);
  template short dot(const svec &v1, const svec &v2);
  template bin dot(const bvec &v1, const bvec &v2);

#if !defined(HAVE_BLAS)
  template double operator*(const vec &v1, const vec &v2);
#if defined(HAVE_MKL) && defined(_MSC_VER)
  template std::complex<double> operator*(const cvec &v1, const cvec &v2);
#endif
#endif
  template int operator*(const ivec &v1, const ivec &v2);
  template short operator*(const svec &v1, const svec &v2);
  template bin operator*(const bvec &v1, const bvec &v2);

#if !defined(HAVE_BLAS)
  template const mat outer_product(const vec &v1, const vec &v2,
				   bool hermitian);
#endif
  template const imat outer_product(const ivec &v1, const ivec &v2,
				    bool hermitian);
  template const smat outer_product(const svec &v1, const svec &v2,
				    bool hermitian);
  template const bmat outer_product(const bvec &v1, const bvec &v2,
				    bool hermitian);

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

  template const vec elem_mult(const vec &a, const vec &b);
  template const cvec elem_mult(const cvec &a, const cvec &b);
  template const ivec elem_mult(const ivec &a, const ivec &b);
  template const svec elem_mult(const svec &a, const svec &b);
  template const bvec elem_mult(const bvec &a, const bvec &b);

  template void elem_mult_out(const vec &a, const vec &b, vec &out);
  template void elem_mult_out(const cvec &a, const cvec &b, cvec &out);
  template void elem_mult_out(const ivec &a, const ivec &b, ivec &out);
  template void elem_mult_out(const svec &a, const svec &b, svec &out);
  template void elem_mult_out(const bvec &a, const bvec &b, bvec &out);

  //------------- Elementwise Multiplication operator (three vectors) ----------

  template const vec elem_mult(const vec &a, const vec &b, const vec &c);
  template const cvec elem_mult(const cvec &a, const cvec &b, const cvec &c);
  template const ivec elem_mult(const ivec &a, const ivec &b, const ivec &c);
  template const svec elem_mult(const svec &a, const svec &b, const svec &c);
  template const bvec elem_mult(const bvec &a, const bvec &b, const bvec &c);

  template void elem_mult_out(const vec &a, const vec &b, const vec &c, vec &out);
  template void elem_mult_out(const cvec &a, const cvec &b, const cvec &c, cvec &out);
  template void elem_mult_out(const ivec &a, const ivec &b, const ivec &c, ivec &out);
  template void elem_mult_out(const svec &a, const svec &b, const svec &c, svec &out);
  template void elem_mult_out(const bvec &a, const bvec &b, const bvec &c, bvec &out);

  //------------- Elementwise Multiplication operator (four vectors) ----------

  template const vec elem_mult(const vec &a, const vec &b, const vec &c, const vec &d);
  template const cvec elem_mult(const cvec &a, const cvec &b, const cvec &c, const cvec &d);
  template const ivec elem_mult(const ivec &a, const ivec &b, const ivec &c, const ivec &d);
  template const svec elem_mult(const svec &a, const svec &b, const svec &c, const svec &d);
  template const bvec elem_mult(const bvec &a, const bvec &b, const bvec &c, const bvec &d);

  template void elem_mult_out(const vec &a, const vec &b, const vec &c, const vec &d, vec &out);
  template void elem_mult_out(const cvec &a, const cvec &b, const cvec &c, const cvec &d, cvec &out);
  template void elem_mult_out(const ivec &a, const ivec &b, const ivec &c, const ivec &d, ivec &out);
  template void elem_mult_out(const svec &a, const svec &b, const svec &c, const svec &d, svec &out);
  template void elem_mult_out(const bvec &a, const bvec &b, const bvec &c, const bvec &d, bvec &out);

  //------------- In-place element-wise multiplication  ----------

  template void elem_mult_inplace(const vec &a, vec &b);
  template void elem_mult_inplace(const cvec &a, cvec &b);
  template void elem_mult_inplace(const ivec &a, ivec &b);
  template void elem_mult_inplace(const svec &a, svec &b);
  template void elem_mult_inplace(const bvec &a, bvec &b);

  //------------- Element-wise multiplication followed by summation ----------

  template double elem_mult_sum(const vec &a, const vec &b);
  template std::complex<double> elem_mult_sum(const cvec &a, const cvec &b);
  template int elem_mult_sum(const ivec &a, const ivec &b);
  template short elem_mult_sum(const svec &a, const svec &b);
  template bin elem_mult_sum(const bvec &a, const bvec &b);

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

  template const vec elem_div(const vec &a, const vec &b);
  template const cvec elem_div(const cvec &a, const cvec &b);
  template const ivec elem_div(const ivec &a, const ivec &b);
  template const svec elem_div(const svec &a, const svec &b);
  template const bvec elem_div(const bvec &a, const bvec &b);

  template const vec elem_div(const double t, const vec &v);
  template const cvec elem_div(const std::complex<double> t, const cvec &v);
  template const ivec elem_div(const int t, const ivec &v);
  template const svec elem_div(const short t, const svec &v);
  template const bvec elem_div(const bin t, const bvec &v);

  template void elem_div_out(const vec &a, const vec &b, vec &out);
  template void elem_div_out(const cvec &a, const cvec &b, cvec &out);
  template void elem_div_out(const ivec &a, const ivec &b, ivec &out);
  template void elem_div_out(const svec &a, const svec &b, svec &out);
  template void elem_div_out(const bvec &a, const bvec &b, bvec &out);

  //------------- Elementwise division followed by summation ----------

  template double elem_div_sum(const vec &a, const vec &b);
  template std::complex<double> elem_div_sum(const cvec &a, const cvec &b);
  template int elem_div_sum(const ivec &a, const ivec &b);
  template short elem_div_sum(const svec &a, const svec &b);
  template bin elem_div_sum(const bvec &a, const bvec &b);

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
