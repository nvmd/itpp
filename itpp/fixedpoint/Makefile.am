include $(top_srcdir)/Makefile.common

h_sources = \
	cfixed.h \
	cfix.h \
	fix_base.h \
	fixed.h \
	fix_factory.h \
	fix_functions.h \
	fix.h \
	fix_operators.h

cpp_sources = \
	cfix.cpp \
	cfixed.cpp \
	fix_base.cpp \
	fix.cpp \
	fixed.cpp \
	fix_factory.cpp \
	fix_functions.cpp \
	fix_operators.cpp

noinst_LTLIBRARIES = libfixedpoint.la
libfixedpoint_la_SOURCES = $(h_sources) $(cpp_sources) 

pkgincludedir = $(includedir)/@PACKAGE@/fixedpoint
pkginclude_HEADERS = $(h_sources)
