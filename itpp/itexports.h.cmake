#ifndef ITEXPORTS_H
#define ITEXPORTS_H

/*define when using the static version of IT++ library*/
#cmakedefine ITPP_STATIC_LIBRARY

#ifndef ITPP_STATIC_LIBRARY

/*needed to export shared library symbols on Windows*/
#ifdef _MSC_VER
  #ifndef ITPP_EXPORT
    #ifdef itpp_EXPORTS /*automatically defined by cmake*/
      #define ITPP_EXPORT __declspec(dllexport)
    #else
      #define ITPP_EXPORT __declspec(dllimport)
    #endif
  #endif
#elif (__GNUC__ >= 4) /*UNIX*/
  #ifndef ITPP_EXPORT
    #define ITPP_EXPORT __attribute__((visibility("default")))
  #endif
#endif

#endif /*ITPP_STATIC_LIBRARY*/

#ifndef ITPP_EXPORT
  #define ITPP_EXPORT
#endif

#endif
