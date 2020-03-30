ACLOCAL_AMFLAGS= -I m4

AM_CFLAGS =  @SIMD_FLAGS@ @OPENMP_CFLAGS@  @TLDEVEL_CPPFLAGS@

CLEANFILES =

AM_LDFLAGS = @DEAD_STRIP@
AM_LDFLAGS += -Wno-undef
AM_LDFLAGS += -static

LIBS = @TLDEVEL_LIB@ $(HDF5_LDFLAGS)  $(HDF5_LIBS) -lm

bin_PROGRAMS = tagcook simdust

simdust_SOURCES = \
sim_dust_test.c \
sim_seq_lib.h \
sim_seq_lib.c

tagcook_SOURCES = \
main.c \
interface.h \
interface.c \
arch_lib.h \
arch_lib.c \
arch_lib_sim.h \
arch_lib_sim.c \
arch_lib_parse_token.h \
arch_lib_parse_token.c \
seq_stats.h \
seq_stats.c \
hmm.h \
hmm.c \
init_hmm.h \
init_hmm.c \
hmm_model_bag.h \
hmm_model_bag.c \
core_hmm_functions.h \
core_hmm_functions.c \
test_arch.h \
test_arch.c \
calibrate_hmm.h \
calibrate_hmm.c \
extract_reads.h \
extract_reads.c \
assign_data.h \
assign_data.c \
assign_actions.c \
filter.h \
filter.c \
bpm.h \
bpm.c \
pst_structs.h \
pst.h \
pst.c \
calibrate_pst.h \
calibrate_pst.c \
lpst.h \
lpst.c \
init_poahmm.h \
poahmm_structs.h \
init_poahmm.c \
poahmm.h \
poahmm.c \
correct.h \
correct.c \
base_quality.h \
base_quality.c \
read_groups.h \
read_groups.c



tagcook_LDADD =  -lpthread -lm

TESTLIST = arch_lib_test interface_test io_hdf5_test zwrite core_hmm_test band_test bq_test read_group_test poahmm_test pst_test correct_test

TESTS = $(TESTLIST)

TESTS_ENVIRONMENT = \
	devdir='$(top_srcdir)'/dev; \
         export devdir; \
	$(VALGRIND)

check_PROGRAMS = $(TESTLIST)

band_test_SOURCES = banded_test.c

core_hmm_test_SOURCES = \
misc.h \
misc.c \
pst_structs.h \
pst.h \
pst.c \
correct.h \
correct.c \
arch_lib.h \
arch_lib.c \
arch_lib_parse_token.h \
arch_lib_parse_token.c \
seq_stats.h \
seq_stats.c \
hmm.h \
init_hmm.h \
init_hmm.c \
hmm_model_bag.h \
hmm_model_bag.c \
core_hmm_functions.h \
core_hmm_functions.c


core_hmm_test_CPPFLAGS = $(AM_CPPFLAGS) -DCORE_HMM_TEST

io_hdf5_test_SOURCES = \
io_hdf5.h \
io_hdf5.c
io_hdf5_test_CPPFLAGS = $(AM_CPPFLAGS)


arch_lib_test_SOURCES = \
pst_structs.h \
pst.h \
pst.c \
correct.h \
correct.c \
arch_lib.h \
arch_lib.c \
arch_lib_parse_token.h \
arch_lib_parse_token.c

arch_lib_test_CPPFLAGS = $(AM_CPPFLAGS) -DARCH_TEST
arch_lib_test_LDADD =

interface_test_SOURCES = \
interface.h \
interface.c
interface_test_CPPFLAGS = $(AM_CPPFLAGS) -DINTERFACE_TEST
interface_test_LDADD =

bq_test_SOURCES = \
base_quality.h \
base_quality.c
bq_test_CPPFLAGS = $(AM_CPPFLAGS) -DBASEQTEST

read_group_test_SOURCES = \
read_groups.h \
read_groups.c \
seq_stats.h \
seq_stats.c 

read_group_test_CPPFLAGS = $(AM_CPPFLAGS) -DREAD_GROUP_ITEST

poahmm_test_SOURCES = \
poahmm_test.c \
init_poahmm.h \
init_poahmm.c \
poahmm.h \
poahmm_structs.h \
poahmm.c \
base_quality.h \
base_quality.c \
sim_seq_lib.h \
sim_seq_lib.c \
arch_lib.h \
arch_lib.c \
pst.h \
pst.c \
correct.h \
correct.c \
arch_lib_parse_token.h \
arch_lib_parse_token.c 


pst_test_SOURCES = \
pst.h \
pst.c \
pst_test.c \
sim_seq_lib.h \
sim_seq_lib.c

correct_test_SOURCES = \
correct.h \
correct.c \
base_quality.h \
base_quality.c \
correct_test.c


zwrite_SOURCES = test_pipe.c

