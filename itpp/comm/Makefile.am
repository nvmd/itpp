include $(top_srcdir)/Makefile.common

h_sources = \
	bch.h \
	channel_code.h \
	channel.h \
	commfunc.h \
	convcode.h \
	crc.h \
	egolay.h \
	error_counters.h \
	galois.h \
	hammcode.h \
	interleave.h \
	modulator.h \
	ofdm.h \
	pulse_shape.h \
	punct_convcode.h \
	rec_syst_conv_code.h \
	reedsolomon.h \
	sequence.h \
	spread.h \
	turbo.h

cpp_sources = \
	bch.cpp \
	channel.cpp \
	commfunc.cpp \
	convcode.cpp \
	crc.cpp \
	egolay.cpp \
	error_counters.cpp \
	galois.cpp \
	hammcode.cpp \
	interleave.cpp \
	modulator.cpp \
	ofdm.cpp \
	pulse_shape.cpp \
	punct_convcode.cpp \
	rec_syst_conv_code.cpp \
	reedsolomon.cpp \
	sequence.cpp \
	spread.cpp \
	turbo.cpp

noinst_LTLIBRARIES = libcomm.la
libcomm_la_SOURCES = $(h_sources) $(cpp_sources) 

pkgincludedir = $(includedir)/@PACKAGE@/comm
pkginclude_HEADERS = $(h_sources)
