include $(top_srcdir)/Makefile.common

SUBDIRS = base comm fixedpoint protocol srccode

h_sources = \
	alpha_modules.h \
	config.h \
	itbase.h \
	itcomm.h \
	itconfig.h \
	itfixedpoint.h \
	itmex.h \
	itprotocol.h \
	itsrccode.h

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
