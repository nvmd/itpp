include $(top_srcdir)/Makefile.common

h_sources = \
	libharu.h \
	plot_pdf.h

cpp_sources = \
	plot_pdf.h

noinst_LTLIBRARIES = libgraphics.la
libgraphics_la_SOURCES = $(h_sources) $(cpp_sources) 

pkgincludedir = $(includedir)/@PACKAGE@/graphics
pkginclude_HEADERS = $(h_sources)
