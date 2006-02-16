include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

SUBDIRS = bessel

noinst_LTLIBRARIES = libbase.la
if ENABLE_DEBUG
noinst_LTLIBRARIES += libbase_debug.la
endif

libbase_la_SOURCES = $(h_base_sources) $(cpp_base_sources)
libbase_la_CXXFLAGS = $(CXXFLAGS_OPT)
libbase_la_LIBADD = $(top_builddir)/@PACKAGE@/base/bessel/libbessel.la

libbase_debug_la_SOURCES = $(h_base_sources) $(cpp_base_sources)
libbase_debug_la_CXXFLAGS = $(CXXFLAGS_DEBUG)
libbase_debug_la_LIBADD = \
	$(top_builddir)/@PACKAGE@/base/bessel/libbessel_debug.la

pkgincludedir = $(includedir)/@PACKAGE@/base
pkginclude_HEADERS = $(h_base_sources)
