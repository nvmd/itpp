include $(top_srcdir)/Makefile.common

ACLOCAL_AMFLAGS = -I config

if HTML_DOCS
HTMLDOCDIR = doc
endif

SUBDIRS = itpp $(HTMLDOCDIR)

EXTRA_DIST = autogen.sh Makefile.common TODO

bin_SCRIPTS = @PACKAGE@-config

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = @PACKAGE@.pc

uninstall-local:
	rm -rf $(DESTDIR)$(includedir)/@PACKAGE@
