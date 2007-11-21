##### http://autoconf-archive.cryp.to/ac_cxx_extern_template.html
#
# SYNOPSIS
#
#   AC_CXX_EXTERN_TEMPLATE
#
# DESCRIPTION
#
#   Test whether the C++ compiler supports "extern template".
#
# LAST MODIFICATION
#
#   2007-11-21  by Adam Piatyszek <ediap@users.sourceforge.net>
#
# COPYLEFT
#
#   Copyright (c) 2005 Patrick Mauritz <oxygene@studentenbude.ath.cx>
#
#   Copying and distribution of this file, with or without
#   modification, are permitted in any medium without royalty provided
#   the copyright notice and this notice are preserved.

AC_DEFUN([AC_CXX_EXTERN_TEMPLATE],
[AC_CACHE_CHECK([whether the compiler supports extern template],
                [ac_cv_cxx_extern_template],
                [AC_LANG_PUSH([C++])
                 AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
[[template <typename T> void foo(T);
extern template void foo<int>(int);]],
                                                    [])],
                                   [ac_cv_cxx_extern_template=yes],
                                   [ac_cv_cxx_extern_template=no])
                 AC_LANG_POP([C++])])
 if test "$ac_cv_cxx_extern_template" = yes; then
   AC_DEFINE([HAVE_EXTERN_TEMPLATE], [1],
             [Define if the compiler supports extern template])
 fi
])
