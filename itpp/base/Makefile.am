include $(top_srcdir)/Makefile.common

SUBDIRS = bessel

h_sources = \
	array.h \
	bessel.h \
	binary.h \
	binfile.h \
	cblas.h \
	cholesky.h \
	circular_buffer.h \
	converters.h \
	copy_vector.h \
	det.h \
	eigen.h \
	elmatfunc.h \
	factory.h \
	fastica.h \
	fastmath.h \
	fftw3.h \
	filter.h \
	freq_filt.h \
	gf2mat.h \
	help_functions.h \
	integration.h \
	inv.h \
	itassert.h \
	itfile.h \
	itpp_version.h \
	lapack.h \
	ls_solve.h \
	lu.h \
	machdep.h \
	matfunc.h \
	mat.h \
	operators.h \
	parser.h \
	qr.h \
	random.h \
	scalfunc.h \
	sigfun.h \
	smat.h \
	sort.h \
	source.h \
	specmat.h \
	stack.h \
	stat.h \
	svd.h \
	svec.h \
	timing.h \
	transforms.h \
	vec.h

cpp_sources = \
	bessel.cpp \
	binary.cpp \
	binfile.cpp \
	cholesky.cpp \
	converters.cpp \
	det.cpp \
	eigen.cpp \
	elmatfunc.cpp \
	fastica.cpp \
	fastmath.cpp \
	filter.cpp \
	freq_filt.cpp \
	gf2mat.cpp \
	integration.cpp \
	inv.cpp \
	itassert.cpp \
	itfile.cpp \
	itpp_version.cpp \
	ls_solve.cpp \
	lu.cpp \
	mat.cpp \
	matfunc.cpp \
	operators.cpp \
	parser.cpp \
	qr.cpp \
	random.cpp \
	scalfunc.cpp \
	sigfun.cpp \
	smat.cpp \
	source.cpp \
	specmat.cpp \
	stat.cpp \
	svd.cpp \
	svec.cpp \
	timing.cpp \
	transforms.cpp \
	vec.cpp

noinst_LTLIBRARIES = libbase.la
libbase_la_SOURCES = $(h_sources) $(cpp_sources)
libbase_la_LIBADD = $(top_builddir)/@PACKAGE@/base/bessel/libbessel.la

pkgincludedir = $(includedir)/@PACKAGE@/base
pkginclude_HEADERS = $(h_sources)
