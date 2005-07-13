include $(top_srcdir)/Makefile.common

h_sources = \
	audiofile.h \
	gmm.h \
	lpcfunc.h \
	pnm.h \
	vq.h \
	vqtrain.h

cpp_sources = \
	audiofile.cpp \
	gmm.cpp \
	lpcfunc.cpp \
	pnm.cpp \
	vq.cpp \
	vqtrain.cpp

noinst_LTLIBRARIES = libsrccode.la
libsrccode_la_SOURCES = $(h_sources) $(cpp_sources) 

pkgincludedir = $(includedir)/@PACKAGE@/srccode
pkginclude_HEADERS = $(h_sources)
