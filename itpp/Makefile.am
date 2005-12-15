include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

SUBDIRS = base comm fixedpoint protocol srccode

EXTRA_DIST = config_msvc.h

lib_LTLIBRARIES = lib@PACKAGE@.la
lib@PACKAGE@_la_SOURCES = 
lib@PACKAGE@_la_LDFLAGS = -version-info $(GENERIC_LIBRARY_VERSION)
lib@PACKAGE@_la_LIBADD = \
	$(top_builddir)/@PACKAGE@/base/libbase.la \
	$(top_builddir)/@PACKAGE@/comm/libcomm.la \
	$(top_builddir)/@PACKAGE@/fixedpoint/libfixedpoint.la \
	$(top_builddir)/@PACKAGE@/protocol/libprotocol.la \
	$(top_builddir)/@PACKAGE@/srccode/libsrccode.la

pkgincludedir = $(includedir)/@PACKAGE@
nobase_pkginclude_HEADERS = $(h_sources)
