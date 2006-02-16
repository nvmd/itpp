include $(top_srcdir)/Makefile.common
LDADD = $(top_builddir)/itpp/libitpp.la

EXTRA_DIST = fastica_test_data.txt parser_test_data.txt $(ref_sources) 

check_PROGRAMS = \
	array_test \
	bch_test \
	bessel_test \
	cholesky_test \
	circular_buffer_test \
	convcode_test \
	det_test \
	eigen_test \
	fastica_test \
	fastmath_test \
	filter_test \
	filter_design_test \
	fix_test \
	freq_filt_test \
	galois_test \
	gf2mat_test \
	integration_test \
	interleaver_test \
	inv_test \
	llr_test \
	ls_solve_test \
	lu_test \
	matfunc_test \
	mat_test \
	modulator_test \
	modulator_nd_test \
	newton_search_test \
	pam_test \
	parser_test \
	poly_test \
	psk_test \
	pulse_shape_test \
	qam_test \
	qr_test \
	rand_test \
	rec_syst_conv_code_test \
	schur_test \
	sigfun_test \
	source_test \
	sparse_test \
	stat_test \
	svd_test \
	timer_test \
	transforms_test \
	turbo_test \
	vec_test \
	vector_test \
	window_test

array_test_SOURCES = array_test.cpp
bch_test_SOURCES = bch_test.cpp
bessel_test_SOURCES = bessel_test.cpp
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
integration_test_SOURCES = integration_test.cpp
interleaver_test_SOURCES = interleaver_test.cpp
inv_test_SOURCES = inv_test.cpp
llr_test_SOURCES = llr_test.cpp
ls_solve_test_SOURCES = ls_solve_test.cpp
lu_test_SOURCES = lu_test.cpp
matfunc_test_SOURCES = matfunc_test.cpp
mat_test_SOURCES = mat_test.cpp
modulator_test_SOURCES = modulator_test.cpp
modulator_nd_test_SOURCES = modulator_nd_test.cpp
newton_search_test_SOURCES = newton_search_test.cpp
pam_test_SOURCES = pam_test.cpp
parser_test_SOURCES = parser_test.cpp
parser_test_CPPFLAGS = -DPARSER_TEST_FILE=\"$(srcdir)/parser_test_data.txt\" $(AM_CPPFLAGS)
poly_test_SOURCES = poly_test.cpp
psk_test_SOURCES = psk_test.cpp
pulse_shape_test_SOURCES = pulse_shape_test.cpp
qam_test_SOURCES = qam_test.cpp
qr_test_SOURCES = qr_test.cpp
rand_test_SOURCES = rand_test.cpp
rec_syst_conv_code_test_SOURCES = rec_syst_conv_code_test.cpp
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
vector_test_SOURCES = vector_test.cpp
window_test_SOURCES = window_test.cpp

noext_tests = $(patsubst %$(EXEEXT), %, $(check_PROGRAMS))
ref_sources = $(patsubst %$(EXEEXT), %.ref, $(check_PROGRAMS))

check-local: $(check_PROGRAMS) $(ref_sources)
	@for i in $(noext_tests); do \
		echo ------------------------------------------------------------------------------; \
		./$$i$(EXEEXT) > $$i.tmp; \
		diff $$i.tmp $(srcdir)/$$i.ref > /dev/null; \
		if test "x$$?" = "x0"; then \
			echo "Test \`$$i' PASSED."; \
		else \
			echo "Test \`$$i' FAILED!!!"; \
		fi; \
	done; \
	echo ------------------------------------------------------------------------------

clean-local:
	rm -f *.tmp
