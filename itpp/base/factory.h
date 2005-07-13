/*---------------------------------------------------------------------------*
 *                                   IT++			             *
 *---------------------------------------------------------------------------*
 * Copyright (c) 2004-2005 by Johan Bergman.                                 *
 *                                                                           *
 * Permission to use, copy, modify, and distribute this software and its     *
 * documentation under the terms of the GNU General Public License is hereby *
 * granted. No representations are made about the suitability of this        *
 * software for any purpose. It is provided "as is" without expressed or     *
 * implied warranty. See the GNU General Public License for more details.    *
 *---------------------------------------------------------------------------*/

/*!
  \file
  \brief Base class for class factories
  \author Johan Bergman

  $Revision$

  $Date$
*/

#ifndef __factory_h
#define __factory_h

#include <cstring>

namespace itpp {

  // Forward declarations
  template<class T> class Array;
  template<class Num_T> class Mat;
  template<class Num_T> class Vec;

  /*!
    \brief Base class for class factories

    A class factory (or virtual constructor) is a class that can create instances
    of another class. Factory is a base class for such factories. When declaring
    an Array, Vec or Mat, a factory can be passed as an (optional) constructor
    argument:
    \code
    // Declare a Vec<type> with size=10 and factory=DEFAULT_FACTORY
    Vec<type> a(10);

    // Declare a Vec<type> with size=10 and factory=f
    Factory f;
    Vec<type> b(10, f);
    \endcode

    By default, the factory (\c DEFAULT_FACTORY and \c f in the above examples)
    is not used at all! However, by overloading a help function called
    \e create_elements we can force Array/Vec/Mat to use the factory for element
    creation (instead of using the default constructor for the element type).

    \note It is the \e numeric elements that will be created by the factory,
    i.e. for an Array<Mat<T> >, the factory will be used for creating the Mat
    elements rather than the Array elements.

    Here is an example that (partly) defines a user-defined numeric type My_Type,
    a corresponding factory My_Factory and a corresponding help function
    create_elements<My_Type> that will be used by Array, Vec and Mat for element
    creation.
    \code
    class My_Type {
    public:
    // Default constructor
    My_Type() : data(0) {}
    // Constructor
    My_Type(int d) : data(d) {}
    .
    .
    .
    protected:
    int data;
    };

    class My_Factory : public Factory {
    public:
    // Constructor
    explicit My_Factory(int d) : init_data(d) {}
    // Destructor
    virtual ~My_Factory() {}
    // Create an n-length array of My_Type
    virtual void create(My_Type* &ptr, const int n) const {ptr = new My_Type[n](init_data);}
    protected:
    int init_data;
    };

    // Create an n-length array of My_Type using My_Factory f
    template<>
    void create_elements<My_Type>(My_Type* &ptr, const int n, const Factory &f)
    {
    if (const My_Factory *my_factory_ptr = dynamic_cast<const My_Factory*>(&f)) {
    // Yes, f seems to be a My_Factory. Now call the My_Factory::create method
    my_factory_ptr->create(ptr, n);
    }
    else {
    // No, f does not seem to be a My_Factory. As a fallback solution,
    // assume that f is DEFAULT_FACTORY and use the default constructor
    ptr = new My_Type[n];
    }
    }
    \endcode

    Now,  
    \code
    // Declare a My_Factory for init_data = 123
    My_Factory my123_factory(123);

    // Declare a Vec<My_Type> with size 10 that uses My_Type() for element creation
    Vec<My_Type> v1(10);

    // Declare a Vec<My_Type> with size 10 that uses My_Type(123) for element creation
    Vec<My_Type> v1(10, my123_factory);
    \endcode

    For a more interesting example, see Fix_Factory.
  */
  class Factory {
  public:
    //! Destructor
    virtual ~Factory() {}
  };

  //! Default (dummy) factory
  const Factory DEFAULT_FACTORY;

  //! Create an n-length array of T to be used as Array, Vec or Mat elements
  template<class T>
  void create_elements(T* &ptr, const int n, const Factory &)
  {
    ptr = new T[n];
  }

  //! Create an n-length array of Array<T> to be used as Array elements
  template<class T>
  void create_elements(Array<T>* &ptr, const int n, const Factory &f)
  {
    ptr = new Array<T>[n];
    if (ptr) {
      // Set factory
      for (int i=0; i<n; i++) {
        // Note: explicit destructor call intentionally omitted (to save time)
        new (ptr + i) Array<T>(f);
      }
    }
  }

  //! Create an n-length array of Mat<T> to be used as Array elements
  template<class T>
  void create_elements(Mat<T>* &ptr, const int n, const Factory &f)
  {
    ptr = new Mat<T>[n];
    if (ptr) {
      // Set factory
      for (int i=0; i<n; i++) {
        // Note: explicit destructor call intentionally omitted (to save time)
        new (ptr + i) Mat<T>(f);
      }
    }
  }

  //! Create an n-length array of Vec<T> to be used as Array elements
  template<class T>
  void create_elements(Vec<T>* &ptr, const int n, const Factory &f)
  {
    ptr = new Vec<T>[n];
    if (ptr) {
      // Set factory
      for (int i=0; i<n; i++) {
        // Note: explicit destructor call intentionally omitted (to save time)
        new (ptr + i) Vec<T>(f);
      }
    }
  }

} //namespace itpp

#endif // __factory_h
