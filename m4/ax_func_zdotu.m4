AC_DEFUN([AX_FUNC_ZDOTU_RETURN],
[AC_CACHE_CHECK([whether zdotu_ returns result],
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
  const int size = 6;
  cdouble a[size] = {cdouble( 0.7219, 0.8871), cdouble( 0.7073,-0.7953),
                     cdouble( 0.2610, 0.4325), cdouble(-0.0565,-0.0719),
                     cdouble( 0.7277,-0.9754), cdouble(-0.3780, 1.0718)};
  cdouble b[size] = {cdouble(-0.0821,+0.8410), cdouble(-0.0749, 0.0729),
                     cdouble(-0.6094,-0.2975), cdouble( 0.2106,-0.2026),
                     cdouble( 0.1043,-0.8300), cdouble( 0.0806, 0.3698)};
  cdouble x_ref(-2.01767031,-0.45861365);
  cdouble x = zdotu_(&size, a, &incr, b, &incr);
  return (std::abs(x - x_ref) < 1e-6) ? 0 : 1;
}]])],
  [ax_cv_zdotu_ret_complex=yes], [ax_cv_zdotu_ret_complex=no])
AC_LANG_POP([C++])])
if test "$ax_cv_zdotu_ret_complex" = yes; then
  AC_DEFINE([HAVE_ZDOTU_RETURN], 1,
            [Define if zdotu_ returns result.])
fi]) # AX_FUNC_ZDOTU_RETURN

AC_DEFUN([AX_FUNC_ZDOTU_VOID],
[AC_CACHE_CHECK([whether zdotu_ passes result as its first argument],
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
  const int size = 6;
  cdouble a[size] = {cdouble( 0.7219, 0.8871), cdouble( 0.7073,-0.7953),
                     cdouble( 0.2610, 0.4325), cdouble(-0.0565,-0.0719),
                     cdouble( 0.7277,-0.9754), cdouble(-0.3780, 1.0718)};
  cdouble b[size] = {cdouble(-0.0821,+0.8410), cdouble(-0.0749, 0.0729),
                     cdouble(-0.6094,-0.2975), cdouble( 0.2106,-0.2026),
                     cdouble( 0.1043,-0.8300), cdouble( 0.0806, 0.3698)};
  cdouble x_ref(-2.01767031,-0.45861365);
  cdouble x;
  zdotu_(&x, &size, a, &incr, b, &incr);
  return (std::abs(x - x_ref) < 1e-6) ? 0 : 1;
}]])],
   [ax_cv_zdotu_ret_void=yes], [ax_cv_zdotu_ret_void=no])
 AC_LANG_POP([C++])])
if test "$ax_cv_zdotu_ret_void" = yes; then
  AC_DEFINE([HAVE_ZDOTU_VOID], 1,
            [Define if zdotu_ passes result as its first argument.])
fi]) # AX_FUNC_ZDOTU_VOID
