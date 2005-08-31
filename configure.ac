# Process this file with autoconf to produce a configure script
# Copyright (c) 2005 by Adam Piatyszek <ediap@users.sourceforge.net>
#-------------------------------------------------------------------
AC_PREREQ(2.59)
AC_INIT([IT++], [3.9.0], [ediap@users.sourceforge.net], [itpp])
AC_CONFIG_SRCDIR([itpp/base/vec.cpp])
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE

# Shared library versioning
GENERIC_LIBRARY_VERSION="1:0:0"
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
CXXFLAGS="-g"
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
  [AS_HELP_STRING([--enable-debug], 
    [enable debugging (disabled by default)])])
test "x$enable_debug" != xyes && enable_debug=no

# Check for pedantic switch
AC_ARG_ENABLE(pedantic, 
  [AS_HELP_STRING([--enable-pedantic], 
    [enable additional warning and pedantic flags (disabled by default)])])
test "x$enable_pedantic" != xyes && enable_pedantic=no

# Check for htmldocs switch
AC_ARG_ENABLE(html-docs, 
  [AS_HELP_STRING([--disable-html-docs], 
    [disable creation of HTML documentation (enabled by default)])])
if test "x$enable_html_docs" = xno; then
  enable_html_docs=disabled
else
  enable_html_docs=yes
fi

# Provide --with-docdir
AC_ARG_WITH(docdir,
	[AS_HELP_STRING([--with-docdir=DIR], [set documentation directory to DIR])],
	[DOC_DIR="$withval"], [DOC_DIR=""])
if test "x$DOC_DIR" != x; then
	docdir="$DOC_DIR"
else
	docdir="${datadir}/doc/$PACKAGE-$VERSION"
fi
AC_SUBST(docdir)

# Check for CPU option
AC_ARG_WITH(cpu, 
  [AS_HELP_STRING([--with-cpu=CPU], 
    [enable optimization for a specific CPU: pentium3, pentium4, athlon-xp, sparc])])

# Set debug/optimised flags
if test "x$enable_debug" = xyes; then
  AC_DEFINE([ASSERT_LEVEL], [2], [Define ASSERT_LEVEL as 2])
  CXXFLAGS="-g -DASSERT_LEVEL=2 $CXXFLAGS"
else
  AC_DEFINE([ASSERT_LEVEL], [1], [Define ASSERT_LEVEL as 1])
  AC_DEFINE([NDEBUG], [1], [Define NDEBUG flag])
  CXXFLAGS="-O3 -DNDEBUG -DASSERT_LEVEL=1 $CXXFLAGS"
fi


# Only for GCC compiler 
if test "x$GXX" = xyes; then
  # Add gcc warnings, in debug/maintainer mode only
  if test "x$enable_debug" = xyes; then
    if test "x$enable_pedantic" = xyes; then
      CXXFLAGS="-W -Wcast-qual -Wcast-align -Wpointer-arith -Wlong-long -Wshadow -Wbad-function-cast -Wwrite-strings -Wstrict-prototypes -Wredundant-decls -Wnested-externs -Woverloaded-virtual -pedantic $CXXFLAGS"
    else
      CXXFLAGS="-Wall $CXXFLAGS"
    fi
  fi

  # Processor flags settings
  case $with_cpu in
    pentium3) CXXFLAGS="-march=pentium3 -ffast-math -mfpmath=sse -msse $CXXFLAGS" ;;
    pentium4) CXXFLAGS="-march=pentium4 -ffast-math -mfpmath=sse -msse2 $CXXFLAGS" ;;
    athlon-xp) CXXFLAGS="-march=athlon-xp -ffast-math -mfpmath=sse -msse2 $CXXFLAGS" ;;
    sparc) CXXFLAGS="-mv8 $CXXFLAGS" ;;
    *) ;;
  esac
fi


# Checks for documentation build tools
AC_CHECK_PROG([doxygen_ok], [doxygen], [yes], [no])
AC_CHECK_PROG([latex_ok], [latex], [yes], [no])
AC_CHECK_PROG([dvips_ok], [dvips], [yes], [no])
AC_CHECK_PROG([gs_ok], [gs], [yes], [no])
if test "x$enable_html_docs" != xdisabled; then
  test "x$doxygen_ok" != xyes && enable_html_docs=no
  test "x$latex_ok" != xyes && enable_html_docs=no
  test "x$dvips_ok" != xyes && enable_html_docs=no
  test "x$gs_ok" != xyes && enable_html_docs=no
  AC_PATH_PROG([PERL], [perl], [])
fi
AM_CONDITIONAL([HTML_DOCS], [test "x$enable_html_docs" = xyes])


# Check for diff program used for tests
AC_CHECK_PROG([diff_ok], [diff], [yes], [no])
AM_CONDITIONAL([HAVE_DIFF], [test "x$diff_ok" = xyes])


# Check for standard math library
AC_CHECK_LIB(m, sin)


# Check F77 names of BLAS and LAPACK functions
AC_F77_FUNC(sgemm) 
AC_F77_FUNC(cheev)

# Checks for Intel MKL
AC_ARG_WITH(mkl_dir, [AS_HELP_STRING([--with-mkl-dir=DIR], 
  [set Intel MKL directory to DIR])], [MKL_DIR=$withval], [MKL_DIR=""])
if test "x$MKL_DIR" != x; then
  MKL_LDFLAGS="-L$MKL_DIR -L$MKL_DIR/lib/32"
fi

mkl_ok=no
AC_ARG_WITH(mkl, [AS_HELP_STRING([--with-mkl=LIB], 
  [use MKL libraries LIB, e.g. -lmkl_p4 -lmkl_lapack -lguide -lpthread])])
case $with_mkl in
  yes | "") ;;
  no) mkl_ok=disabled ;;
  *) MKL_LIBS="$with_mkl" ;;
esac

if test "x$mkl_ok" = xno; then
  if test "x$MKL_LIBS" != x; then
    mkl_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$MKL_LDFLAGS $LDFLAGS"
    mkl_save_LIBS="$LIBS"; LIBS="$MKL_LIBS $LIBS $FLIBS"
    AC_CHECK_FUNC($sgemm,
      [AC_CHECK_FUNC($cheev, 
        [AC_CHECK_FUNC(cblas_sgemm,
          [AC_CHECK_FUNC(DftiComputeForward, [mkl_ok=yes])])])])
    LIBS="$mkl_save_LIBS"
    LDFLAGS="$mkl_save_LDFLAGS"
  fi
fi

if test "x$mkl_ok" = xyes; then
  AC_DEFINE([HAVE_MKL], [1], [Define if you have MKL library.])
  LIBS="$MKL_LIBS $LIBS $FLIBS"
  LDFLAGS="$MKL_LDFLAGS $LDFLAGS"
  blas_ok=disabled
  atlas_ok=disabled
  fftw_ok=disabled
  cblas_ok=disabled
  lapack_ok=disabled
fi


# Checks for ATLAS and BLAS libraries
if test "x$atlas_ok" != xdisabled; then
  atlas_save_LIBS="$LIBS"; LIBS="$LIBS $FLIBS"
  atlas_ok=no
  test "x$blas_ok" != xdisabled && blas_ok=no
  AC_CHECK_LIB(atlas, ATL_xerbla,
    [AC_CHECK_LIB(f77blas, $sgemm,
      [AC_CHECK_LIB(cblas, cblas_dgemm, 
        [atlas_ok=yes; blas_ok=yes; ATLAS_LIBS="-lcblas -lf77blas -latlas"; 
          BLAS_LIBS="-lf77blas"],
        [], [-lf77blas -latlas])],
      [AC_CHECK_LIB(blas, $sgemm,
        [AC_CHECK_LIB(cblas, cblas_dgemm,
          [atlas_ok=yes; blas_ok=yes; ATLAS_LIBS="-lcblas -lblas -latlas";
            BLAS_LIBS="-lblas"],
          [], [-lblas -latlas])], 
        [], [-latlas])])],
    [AC_CHECK_LIB(blas, $sgemm, [blas_ok=yes; BLAS_LIBS="-lblas"])])
  LIBS="$atlas_save_LIBS"
fi

# Checks for CBLAS library
AC_ARG_WITH(cblas_dir, [AS_HELP_STRING([--with-cblas-dir=DIR], 
  [set CBLAS library directory to DIR])], [CBLAS_DIR=$withval], [CBLAS_DIR=""])
if test "x$CBLAS_DIR" != x; then
  CBLAS_LDFLAGS="-L$CBLAS_DIR"
fi
cblas_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$CBLAS_LDFLAGS $LDFLAGS"

test "x$cblas_ok" != xdisabled && cblas_ok=no
AC_ARG_WITH(cblas, [AS_HELP_STRING([--with-cblas=LIB], 
  [use CBLAS library, e.g. -lcblas -lblas])])
case $with_cblas in
  yes | "") ;;
  no) cblas_ok=disabled ;;
  *) CBLAS_LIBS="$with_cblas" ;;
esac

if test "x$atlas_ok" = xyes; then
  cblas_ok=yes
  CBLAS_LIBS="$ATLAS_LIBS"
fi

if test "x$cblas_ok" = xno; then
  if test "x$CBLAS_LIBS" != x; then
    cblas_save_LIBS="$LIBS"; LIBS="$CBLAS_LIBS $LIBS $FLIBS"
    AC_CHECK_FUNC(cblas_sgemm, [cblas_ok=yes], [CBLAS_LIBS=""])
    LIBS="$cblas_save_LIBS"
  fi
fi

if test "x$cblas_ok" = xno; then
  cblas_save_LIBS="$LIBS"; LIBS="$LIBS $FLIBS"
  AC_CHECK_LIB(cblas, cblas_sgemm, [cblas_ok=yes; CBLAS_LIBS="-lcblas -lblas"],
    [AC_CHECK_LIB(blas, cblas_sgemm, [cblas_ok=yes; CBLAS_LIBS="-lblas"])], 
    [-lblas])
  LIBS="$cblas_save_LIBS"
fi

if test "x$cblas_ok" = xyes; then
  AC_DEFINE([HAVE_CBLAS], [1], [Define if you have CBLAS library.])
  LIBS="$CBLAS_LIBS $LIBS $FLIBS"
else
  LDFLAGS="$cblas_save_LDFLAGS"
fi


# Checks for LAPACK library
AC_ARG_WITH(lapack_dir, [AS_HELP_STRING([--with-lapack-dir=DIR], 
  [set LAPACK library directory to DIR])], 
  [LAPACK_DIR=$withval], [LAPACK_DIR=""])
if test "x$LAPACK_DIR" != x; then
  LAPACK_LDFLAGS="-L$LAPACK_DIR"
fi
lapack_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$LAPACK_LDFLAGS $LDFLAGS"

test "x$lapack_ok" != xdisabled && lapack_ok=no
AC_ARG_WITH(lapack, [AS_HELP_STRING([--with-lapack=LIB], 
  [use LAPACK library, e.g. -llapack -lblas])])
case $with_lapack in
  yes | "") ;;
  no) lapack_ok=disabled ;;
  *) LAPACK_LIBS="$with_lapack" ;;
esac

if test "x$lapack_ok" = xno; then
  if test "x$LAPACK_LIBS" != x; then
    lapack_save_LIBS="$LIBS"; LIBS="$LAPACK_LIBS $LIBS $FLIBS"
    AC_CHECK_FUNC($cheev, [lapack_ok=yes], [LAPACK_LIBS=""])
    LIBS="$lapack_save_LIBS"
  fi
fi

if test "x$blas_ok" = xyes; then
  if test "x$atlas_ok" = xyes; then
    BLAS_LIBS="$ATLAS_LIBS"
  fi
  if test "x$lapack_ok" = xno; then
    lapack_save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS $FLIBS"
    AC_CHECK_FUNC($cheev, [lapack_ok=yes])
    LIBS="$lapack_save_LIBS"
  fi
  if test "x$lapack_ok" = xno; then
    for lapack in lapack lapack_rs6k; do
      lapack_save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS $FLIBS"
      AC_CHECK_LIB($lapack, $cheev, [lapack_ok=yes; LAPACK_LIBS="-l$lapack"])
      LIBS="$lapack_save_LIBS"
    done
  fi
fi
    
if test "x$lapack_ok" = xyes; then
  AC_DEFINE([HAVE_LAPACK], [1], [Define if you have LAPACK library.])
  if test "x$cblas_ok" = xyes; then
    LIBS="$LAPACK_LIBS $LIBS"
  else
    LIBS="$LAPACK_LIBS $BLAS_LIBS $LIBS $FLIBS"
  fi
else    
  LDFLAGS="$lapack_save_LDFLAGS"
fi


# Checks for FFTW library
AC_ARG_WITH(fftw_dir, [AS_HELP_STRING([--with-fftw-dir=DIR], 
  [set FFTW library directory to DIR])], [FFTW_DIR=$withval], [FFTW_DIR=""])
if test "x$FFTW_DIR" != x; then
  FFTW_LDFLAGS="-L$FFTW_DIR"
fi
fftw_save_LDFLAGS="$LDFLAGS"; LDFLAGS="$FFTW_LDFLAGS $LDFLAGS"

test "x$fftw_ok" != xdisabled && fftw_ok=no
AC_ARG_WITH(fftw, [AS_HELP_STRING([--with-fftw=LIB], 
  [use FFTW library, e.g. -lfftw3])])
case $with_fftw in
  yes | "") ;;
  no) fftw_ok=disabled ;;
  *) FFTW_LIBS="$with_fftw" ;;
esac

if test "x$fftw_ok" = xno; then
  if test "x$FFTW_LIBS" != x; then
    fftw_save_LIBS="$LIBS"; LIBS="$FFTW_LIBS $LIBS"
    AC_CHECK_FUNC(fftw_plan_dft_1d, [fftw_ok=yes], [FFTW_LIBS=""])
    LIBS="$fftw_save_LIBS"
  fi
fi

if test "x$fftw_ok" = xno; then
  AC_CHECK_LIB(fftw3, fftw_plan_dft_1d, [fftw_ok=yes; FFTW_LIBS="-lfftw3"],
    [], [-lm])
fi

if test "x$fftw_ok" = xyes; then
  AC_DEFINE([HAVE_FFTW], [1], [Define if you have FFTW library.])
  LIBS="$FFTW_LIBS $LIBS"
else
  LDFLAGS="$fftw_save_LDFLAGS"
fi


# Checks for header files
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h string.h sys/time.h unistd.h])

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
  doc/doxygen_html.cfg \
  tests/Makefile \
  itpp.pc \
  itpp-config])

AC_OUTPUT

echo "
------------------------------------------------------------------------------
$PACKAGE-$VERSION library configuration:
------------------------------------------------------------------------------

Directories:
  - prefix ....... : ${prefix}
  - exec_prefix .. : ${exec_prefix}
  - includedir ... : ${includedir}
  - libdir ....... : ${libdir}
  - docdir ....... : ${docdir}

Platform:
  - CPU .......... : ${with_cpu}

Switches:
  - debug ........ : ${enable_debug}
  - pedantic ..... : ${enable_pedantic}
  - html-docs .... : ${enable_html_docs}

Documentation tools:
  - Doxygen ...... : ${doxygen_ok}
  - LaTeX ........ : ${latex_ok}
  - Dvips ........ : ${dvips_ok}
  - GhostScript .. : ${gs_ok}

Testing tools:
  - Diff ......... : ${diff_ok}

External libs:
  - MKL .......... : ${mkl_ok}
  - FFTW ......... : ${fftw_ok}
  - ATLAS ........ : ${atlas_ok}
  - CBLAS ........ : ${cblas_ok}
  - LAPACK ....... : ${lapack_ok}

Compiler/linker flags/libs/defs:
  - CXX .......... : ${CXX}
  - CPPFLAGS ..... : ${CPPFLAGS}
  - CXXFLAGS ..... : ${CXXFLAGS}
  - LDFLAGS ...... : ${LDFLAGS}
  - LIBS ......... : ${LIBS}
  - DEFS ......... : ${DEFS}

------------------------------------------------------------------------------
Now type 'make && make install' to build and install $PACKAGE-$VERSION library
------------------------------------------------------------------------------
"
