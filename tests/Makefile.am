include $(top_srcdir)/Makefile.common
if ENABLE_DEBUG
LDADD = $(top_builddir)/itpp/libitpp_debug.la
AM_CXXFLAGS = $(CXXFLAGS_DEBUG)
else
LDADD = $(top_builddir)/itpp/libitpp.la
AM_CXXFLAGS = $(CXXFLAGS_OPT)
endif

BASE_TESTS = \
	array_test \
	bessel_test \
	cholesky_test \
	circular_buffer_test \
	det_test \
	eigen_test \
	fastmath_test \
	gf2mat_test \
	integration_test \
	inv_test \
	itfile_test \
	ls_solve_test \
	lu_test \
	matfunc_test \
	mat_test \
	parser_test \
	qr_test \
	rand_test \
	schur_test \
	sparse_test \
	svd_test \
	timer_test \
	vec_test

if ENABLE_COMM
  COMM_TESTS = \
	bch_test \
	channel_test \
	convcode_test \
	galois_test \
	interleaver_test \
	ldpc_test \
	llr_test \
	modulator_test \
	modulator_nd_test \
	pulse_shape_test \
	rec_syst_conv_code_test \
	reedsolomon_test \
	turbo_test
endif

if ENABLE_FIXED
  FIXED_TESTS = \
	fix_test
endif

if ENABLE_OPTIM
  OPTIM_TESTS = \
	newton_search_test
endif

if ENABLE_SIGNAL
  SIGNAL_TESTS = \
	fastica_test \
	filter_test \
	filter_design_test \
	freq_filt_test \
	poly_test \
	sigfun_test \
	source_test \
	transforms_test \
	window_test
endif

STAT_TESTS = \
	histogram_test \
	stat_test

ALL_TESTS = $(BASE_TESTS) $(COMM_TESTS) $(FIXED_TESTS) $(OPTIM_TESTS) \
	$(SIGNAL_TESTS) $(STAT_TESTS)


array_test_SOURCES = array_test.cpp
bch_test_SOURCES = bch_test.cpp
bessel_test_SOURCES = bessel_test.cpp
channel_test_SOURCES = channel_test.cpp
itfile_test_SOURCES = itfile_test.cpp
itfile_test_CPPFLAGS = -DITFILE_TEST_FILE=\"$(srcdir)/itfile_test_data.it\" $(AM_CPPFLAGS)
cholesky_test_SOURCES = cholesky_test.cpp
circular_buffer_test_SOURCES = circular_buffer_test.cpp
convcode_test_SOURCES = convcode_test.cpp
det_test_SOURCES = det_test.cpp
eigen_test_SOURCES = eigen_test.cpp
fastica_test_SOURCES = fastica_test.cpp
fastica_test_CPPFLAGS = -DFASTICA_TEST_FILE=\"$(srcdir)/fastica_test_data.txt\" $(AM_CPPFLAGS)
fastmath_test_SOURCES = fastmath_test.cpp
filter_test_SOURCES = filter_test.cpp
filter_design_test_SOURCES = filter_design_test.cpp
fix_test_SOURCES = fix_test.cpp
freq_filt_test_SOURCES = freq_filt_test.cpp
galois_test_SOURCES = galois_test.cpp
gf2mat_test_SOURCES = gf2mat_test.cpp
histogram_test_SOURCES = histogram_test.cpp
integration_test_SOURCES = integration_test.cpp
interleaver_test_SOURCES = interleaver_test.cpp
inv_test_SOURCES = inv_test.cpp
ldpc_test_SOURCES = ldpc_test.cpp
llr_test_SOURCES = llr_test.cpp
ls_solve_test_SOURCES = ls_solve_test.cpp
lu_test_SOURCES = lu_test.cpp
matfunc_test_SOURCES = matfunc_test.cpp
mat_test_SOURCES = mat_test.cpp
modulator_test_SOURCES = modulator_test.cpp
modulator_nd_test_SOURCES = modulator_nd_test.cpp
newton_search_test_SOURCES = newton_search_test.cpp
parser_test_SOURCES = parser_test.cpp
parser_test_CPPFLAGS = -DPARSER_TEST_FILE=\"$(srcdir)/parser_test_data.txt\" $(AM_CPPFLAGS)
poly_test_SOURCES = poly_test.cpp
pulse_shape_test_SOURCES = pulse_shape_test.cpp
qr_test_SOURCES = qr_test.cpp
rand_test_SOURCES = rand_test.cpp
rec_syst_conv_code_test_SOURCES = rec_syst_conv_code_test.cpp
reedsolomon_test_SOURCES = reedsolomon_test.cpp
schur_test_SOURCES = schur_test.cpp
sigfun_test_SOURCES = sigfun_test.cpp
source_test_SOURCES = source_test.cpp
sparse_test_SOURCES = sparse_test.cpp
stat_test_SOURCES = stat_test.cpp
svd_test_SOURCES = svd_test.cpp
timer_test_SOURCES = timer_test.cpp
transforms_test_SOURCES = transforms_test.cpp
turbo_test_SOURCES = turbo_test.cpp
vec_test_SOURCES = vec_test.cpp
window_test_SOURCES = window_test.cpp


check_PROGRAMS = $(ALL_TESTS)
ref_files = $(ALL_TESTS:=.ref)
tmp_files = $(ALL_TESTS:=.tmp)

EXTRA_DIST = itfile_test_data.it fastica_test_data.txt parser_test_data.txt \
	$(ref_files) 

check-local: $(check_PROGRAMS) $(ref_files)
	@for i in $(ALL_TESTS); do \
		echo ------------------------------------------------------------------------------; \
		./$$i$(EXEEXT) 2> /dev/null | sed 's/^M$$//' | sed 's/\(-\?[0-9]\+\.\?[0-9]*e[+-]\)0\([0-9][0-9]\)/\1\2/g' > $$i.tmp; \
		diff $$i.tmp $(srcdir)/$$i.ref > /dev/null; \
		if test "x$$?" = "x0"; then \
			echo "Test \`$$i' PASSED."; \
		else \
			echo "Test \`$$i' FAILED!!!"; \
		fi; \
	done; \
	echo ------------------------------------------------------------------------------

clean-local:
	-rm -f $(tmp_files) gf2mat_test.alist ldpc_test.codec
