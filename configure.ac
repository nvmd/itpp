# Process this file with autoconf to produce a configure script
#-------------------------------------------------------------------
AC_PREREQ(2.59)
AC_COPYRIGHT([Copyright (C) 2005  Adam Piatyszek <ediap@users.sourceforge.net>])
AC_INIT([IT++], [3.10.0], [ediap@users.sourceforge.net], [itpp])
AC_CONFIG_HEADER([itpp/config.h])
AC_CONFIG_SRCDIR([itpp/base/vec.cpp])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE

# RPM release number
RPM_RELEASE="1"
AC_SUBST(RPM_RELEASE)

# Shared library versioning
GENERIC_LIBRARY_VERSION="2:0:0"
#                        | | |
#                 +------+ | +---+
#                 |        |     |
#              current:revision:age
#                 |        |     |
#                 |        |     +- increment if interfaces have been added;
#                 |        |        set to zero if interfaces have been 
#                 |        |        removed or changed
#                 |        +- increment if source code has changed;
#                 |           set to zero if current is incremented
#                 +- increment if interfaces have been added, removed
#                    or changed
AC_SUBST(GENERIC_LIBRARY_VERSION)

# Checks for programs
save_CXXFLAGS="$CXXFLAGS"
AC_PROG_CXX
AC_PROG_CXXCPP
CXXFLAGS="$save_CXXFLAGS"
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_PROG_LN_S
AC_PROG_MAKE_SET

# Set default language to C++
AC_LANG([C++])

# Check for debug switch
AC_ARG_ENABLE(debug,
              [AC_HELP_STRING([--enable-debug],
                              [enable building debug library [default=no]])])
test "x$enable_debug" != xyes && enable_debug=no
AM_CONDITIONAL([ENABLE_DEBUG], [test "x$enable_debug" = xyes])

# Check for exceptions switch
AC_ARG_ENABLE(exceptions,
              [AC_HELP_STRING([--enable-exceptions],
                              [enable exceptions handling [default=no]])],
              [AC_DEFINE([ITPP_EXCEPTIONS], [1],
                         [Define if you want exceptions handling])])
test "x$enable_exceptions" != xyes && enable_exceptions=no

# Check for htmldocs switch
AC_ARG_ENABLE(html-doc, 
              [AC_HELP_STRING([--enable-html-doc],
                              [enable HTML documentation [default=yes]])])
test "x$enable_html_doc" != xno && enable_html_doc=yes

# Provide --with-docdir
AC_ARG_WITH(docdir,
            [AC_HELP_STRING([--with-docdir=DIR],
                            [set documentation directory to DIR])],
            [DOC_DIR="$withval"],
            [DOC_DIR=""])
if test "x$DOC_DIR" != x; then
  docdir="$DOC_DIR"
else
  docdir="${datadir}/doc/$PACKAGE-$VERSION"
fi
AC_SUBST(docdir)


# Set debug flags
if test "x$enable_debug" = xyes; then
  AC_DEFINE([ASSERT_LEVEL], [2], [Define ASSERT_LEVEL as 2])
  CXXFLAGS_DEBUG="-g"
  # Only for GCC compiler 
  if test "x$GXX" = xyes; then
    CXXFLAGS_DEBUG="-Wall $CXXFLAGS_DEBUG"
  fi
else
  AC_DEFINE([ASSERT_LEVEL], [1], [Define ASSERT_LEVEL as 1])
  AC_DEFINE([NDEBUG], [1], [Define NDEBUG flag])
  CXXFLAGS_DEBUG=""
fi

# Set default CXXFLAGS if not set
if test "x$CXXFLAGS" = x; then
  # Only for GCC compiler 
  if test "x$GXX" = xyes; then
    CXXFLAGS_OPT="-O3 -pipe"
  else
    CXXFLAGS_OPT="-O3" 
  fi
else
  CXXFLAGS_OPT="$CXXFLAGS"
  CXXFLAGS=""
fi

AC_SUBST(CXXFLAGS_OPT)
AC_SUBST(CXXFLAGS_DEBUG)


# Checks for documentation build tools
AC_CHECK_PROG([doxygen_ok], [doxygen], [yes], [no])
AC_CHECK_PROG([latex_ok], [latex], [yes], [no])
AC_CHECK_PROG([dvips_ok], [dvips], [yes], [no])
AC_CHECK_PROG([gs_ok], [gs], [yes], [no])
if test "x$enable_html_doc" != xno; then
  test "x$doxygen_ok" != xyes && enable_html_doc=no
  test "x$latex_ok" != xyes && enable_html_doc=no
  test "x$dvips_ok" != xyes && enable_html_doc=no
  test "x$gs_ok" != xyes && enable_html_doc=no
  AC_PATH_PROG([PERL], [perl], [])
fi
AM_CONDITIONAL([HTML_DOC], [test "x$enable_html_doc" = xyes])


# Check for diff program used for tests
AC_CHECK_PROG([diff_ok], [diff], [yes], [no])
AM_CONDITIONAL([HAVE_DIFF], [test "x$diff_ok" = xyes])


# Checks for BLAS (ACML, MKL, ATLAS or NETLIB)
ACX_BLAS

# Checks for LAPACK (ACML, MKL, ATLAS or NETLIB)
ACX_LAPACK

# Checks for CBLAS (MKL, ATLAS or NETLIB)
ACX_CBLAS

# Checks for FFT (MKL or FFTW3)
ACX_FFT


# Set LIBS
if test x"$acx_blas_ok" = xyes; then
  LIBS="$BLAS_LIBS $LIBS $FLIBS"
  if test x"$LAPACK_LIBS" != x; then
    LIBS="$LAPACK_LIBS $LIBS"
  fi
  if test "x$CBLAS_LIBS" != x; then
    LIBS="$CBLAS_LIBS $LIBS"
  fi
fi
if test "x$FFT_LIBS" != x; then
  LIBS="$FFT_LIBS $LIBS"
fi


# Checks for header files
AC_HEADER_STDC
AC_HEADER_STAT
AC_HEADER_TIME
AC_CHECK_HEADERS([cmath \
  complex \
  deque \
  fstream \
  iomanip \
  iostream \
  limits \
  list \
  queue \
  stdexcept \
  string])


# Define configured files
AC_CONFIG_FILES([Makefile \
  itpp/Makefile \
  itpp/base/Makefile \
  itpp/base/bessel/Makefile \
  itpp/comm/Makefile \
  itpp/fixedpoint/Makefile \
  itpp/protocol/Makefile \
  itpp/srccode/Makefile \
  doc/Makefile \
  doc/images/Makefile \
  doc/tutorial/Makefile \
  doc/tutorial/src/Makefile \
  doc/local/Makefile \
  doc/local/index.doc \
  doc/doxygen_html.cfg \
  extras/Makefile \
  tests/Makefile \
  itpp.pc \
  itpp-config \
  itpp.spec])

AC_OUTPUT

echo "
------------------------------------------------------------------------------
$PACKAGE-$VERSION library configuration:
------------------------------------------------------------------------------

Directories:
  - prefix ......... : ${prefix}
  - exec_prefix .... : ${exec_prefix}
  - includedir ..... : ${includedir}
  - libdir ......... : ${libdir}
  - docdir ......... : ${docdir}

Switches:
  - debug .......... : ${enable_debug}
  - exceptions ..... : ${enable_exceptions}
  - html-doc ....... : ${enable_html_doc}
  - shared ......... : ${enable_shared}
  - static ......... : ${enable_static}

Documentation tools:
  - doxygen ........ : ${doxygen_ok}
  - latex .......... : ${latex_ok}
  - dvips .......... : ${dvips_ok}
  - ghostscript .... : ${gs_ok}

Testing tools:
  - diff ........... : ${diff_ok}

External libs:
  - BLAS ........... : ${acx_blas_ok}
    * MKL .......... : ${blas_mkl_ok}
    * ACML ......... : ${blas_acml_ok}
    * ATLAS ........ : ${blas_atlas_ok}
  - CBLAS .......... : ${acx_cblas_ok}
  - LAPACK ......... : ${acx_lapack_ok}
  - FFT ............ : ${acx_fft_ok}
    * MKL .......... : ${fft_mkl8_ok}
    * ACML ......... : ${fft_acml_ok}
    * FFTW ......... : ${fftw3_ok}

Compiler/linker flags/libs/defs:
  - CXX ............ : ${CXX}
  - CXXFLAGS ....... : ${CXXFLAGS_OPT}
  - CXXFLAGS_DEBUG . : ${CXXFLAGS_DEBUG}
  - CPPFLAGS ....... : ${CPPFLAGS}
  - LDFLAGS ........ : ${LDFLAGS}
  - LIBS ........... : ${LIBS}

------------------------------------------------------------------------------
Now type 'make && make install' to build and install $PACKAGE-$VERSION library
------------------------------------------------------------------------------
"
