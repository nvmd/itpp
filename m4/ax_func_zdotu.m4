AC_DEFUN([AX_FUNC_ZDOTU_RETURN],
[AC_CACHE_CHECK([zdotu_ can return result like a function],
                [ax_cv_zdotu_ret_complex],
[AC_LANG_PUSH([C++])
 AC_RUN_IFELSE([AC_LANG_SOURCE(
[[#include <complex>
typedef std::complex<double> cdouble;
extern "C" {
  cdouble zdotu_(const int *, const cdouble *, const int *,
                 const cdouble *, const int *);
}
int main() {
  const int incr = 1;
  const int size = 3;
  cdouble a[size], b[size];
  for (int i = 0; i < size; ++i) {
    a[i] = cdouble(i+1.5, -i);
    b[i] = cdouble(2, 0.5-i);
  }
  cdouble x_ref(11.5, -11.75);
  cdouble x = zdotu_(&size, a, &incr, b, &incr);
  return (x == x_ref) ? 0 : 1;
}]])],
  [ax_cv_zdotu_ret_complex=yes], [ax_cv_zdotu_ret_complex=no])
AC_LANG_POP([C++])])
if test "$ax_cv_zdotu_ret_complex" = yes; then
  AC_DEFINE([HAVE_ZDOTU_RETURN], 1,
            [Define if zdotu_ can return result like a function.])
fi]) # AX_FUNC_ZDOTU_RETURN

AC_DEFUN([AX_FUNC_ZDOTU_VOID],
[AC_CACHE_CHECK([zdotu_ can pass result as its first argument],
                [ax_cv_zdotu_ret_void],
[AC_LANG_PUSH([C++])
 AC_RUN_IFELSE([AC_LANG_SOURCE(
[[#include <complex>
typedef std::complex<double> cdouble;
extern "C" {
  void zdotu_(cdouble *, const int *, const cdouble *, const int *,
              const cdouble *, const int *);
}
int main() {
  const int incr = 1;
  const int size = 3;
  cdouble a[size], b[size];
  for (int i = 0; i < size; ++i) {
    a[i] = cdouble(i+1.5, -i);
    b[i] = cdouble(2, 0.5-i);
  }
  cdouble x_ref(11.5, -11.75);
  cdouble x;
  zdotu_(&x, &size, a, &incr, b, &incr);
  return (x == x_ref) ? 0 : 1;
}]])],
   [ax_cv_zdotu_ret_void=yes], [ax_cv_zdotu_ret_void=no])
 AC_LANG_POP([C++])])
if test "$ax_cv_zdotu_ret_void" = yes; then
  AC_DEFINE([HAVE_ZDOTU_VOID], 1,
            [Define if zdotu_ can pass result as its first argument.])
fi]) # AX_FUNC_ZDOTU_VOID
