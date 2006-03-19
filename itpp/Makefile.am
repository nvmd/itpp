include $(top_srcdir)/Makefile.common
include $(srcdir)/sources.mk

SUBDIRS = base comm fixed protocol srccode

EXTRA_DIST = config_msvc.h

lib_LTLIBRARIES = lib@PACKAGE@.la
if ENABLE_DEBUG
lib_LTLIBRARIES += lib@PACKAGE@_debug.la
endif

lib@PACKAGE@_la_SOURCES = 
lib@PACKAGE@_la_LDFLAGS = -version-info $(GENERIC_LIBRARY_VERSION)
lib@PACKAGE@_la_LIBADD = \
	$(top_builddir)/@PACKAGE@/base/libbase.la \
	$(top_builddir)/@PACKAGE@/comm/libcomm.la \
	$(top_builddir)/@PACKAGE@/fixed/libfixed.la \
	$(top_builddir)/@PACKAGE@/protocol/libprotocol.la \
	$(top_builddir)/@PACKAGE@/srccode/libsrccode.la

lib@PACKAGE@_debug_la_SOURCES = 
lib@PACKAGE@_debug_la_LDFLAGS = -version-info $(GENERIC_LIBRARY_VERSION)
lib@PACKAGE@_debug_la_LIBADD = \
	$(top_builddir)/@PACKAGE@/base/libbase_debug.la \
	$(top_builddir)/@PACKAGE@/comm/libcomm_debug.la \
	$(top_builddir)/@PACKAGE@/fixed/libfixed_debug.la \
	$(top_builddir)/@PACKAGE@/protocol/libprotocol_debug.la \
	$(top_builddir)/@PACKAGE@/srccode/libsrccode_debug.la

pkgincludedir = $(includedir)/@PACKAGE@
nobase_pkginclude_HEADERS = $(h_sources)
