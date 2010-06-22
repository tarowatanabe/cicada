
AM_CPPFLAGS = -I.. $(ICU_CPPFLAGS) $(BOOST_CPPFLAGS) $(PTHREAD_FLAGS) $(TCMALLOC_CPPFLAGS) $(PROFILER_CPPFLAGS) $(SPARSEHASH_CPPFLAGS)

perftools_LDADD = $(JEMALLOC_LDFLAGS) $(TCMALLOC_LDFLAGS) $(PROFILER_LDFLAGS)

boost_LDADD = \
	$(BOOST_THREAD_LDFLAGS) $(BOOST_THREAD_LIBS) \
	$(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
	$(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_IOSTREAMS_LIBS) \
	$(BOOST_SYSTEM_LDFLAGS) $(BOOST_SYSTEM_LIBS)

LIBCICADA = $(top_builddir)/cicada/libcicada.la
LIBUTILS   = $(top_builddir)/utils/libutils.la

bin_PROGRAMS = \
	cicada \
	cicada_mpi \
	cicada_eval \
	cicada_learn \
	cicada_learn_mpi \
	cicada_mert \
	cicada_mert_mpi \
	cicada_maxlike \
	filter_penntreebank \
	filter_cabocha \
	grammar_index \
	mpish \
	mpipe

cicada_SOURCES = cicada.cpp
cicada_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mpi_SOURCES  = cicada_mpi.cpp
cicada_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_eval_SOURCES = cicada_eval.cpp
cicada_eval_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_SOURCES = cicada_learn.cpp lbfgs.c
cicada_learn_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_mpi_SOURCES  = cicada_learn_mpi.cpp lbfgs.c
cicada_learn_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_learn_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mert_SOURCES = cicada_mert.cpp
cicada_mert_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mert_mpi_SOURCES  = cicada_mert_mpi.cpp
cicada_mert_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_mert_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_maxlike_SOURCES = cicada_maxlike.cpp lbfgs.c
cicada_maxlike_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

grammar_index_SOURCES = grammar_index.cpp
grammar_index_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

filter_penntreebank_SOURCES = filter_penntreebank.cpp
filter_penntreebank_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

filter_cabocha_SOURCES = filter_cabocha.cpp
filter_cabocha_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)


mpish_SOURCES = mpish.cpp
mpish_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
mpish_LDADD   = $(MPI_LDFLAGS) $(boost_LDADD) $(perftools_LDADD)

mpipe_SOURCES = mpipe.cpp
mpipe_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
mpipe_LDADD   = $(MPI_LDFLAGS) $(boost_LDADD) $(perftools_LDADD)
