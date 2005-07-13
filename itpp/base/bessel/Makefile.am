include $(top_srcdir)/Makefile.common

h_sources = \
	bessel_internal.h

cpp_sources = \
	airy.cpp \
	chbevl.cpp \
	hyperg.cpp \
	i0.cpp \
	i1.cpp \
	iv.cpp \
	jv.cpp \
	k0.cpp \
	k1.cpp \
	kn.cpp \
	polevl.cpp \
	struve.cpp

noinst_LTLIBRARIES = libbessel.la
libbessel_la_SOURCES = $(h_sources) $(cpp_sources)

pkgincludedir = $(includedir)/@PACKAGE@/base/bessel
pkginclude_HEADERS = $(h_sources)
