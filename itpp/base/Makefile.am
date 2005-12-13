include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

SUBDIRS = bessel

noinst_LTLIBRARIES = libbase.la
libbase_la_SOURCES = $(h_base_sources) $(cpp_base_sources)
libbase_la_LIBADD = $(top_builddir)/@PACKAGE@/base/bessel/libbessel.la

pkgincludedir = $(includedir)/@PACKAGE@/base
pkginclude_HEADERS = $(h_base_sources)
