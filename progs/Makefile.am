
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
	cicada_extract_ghkm \
	cicada_extract_ghkm_mpi \
	cicada_extract_phrase \
	cicada_extract_phrase_mpi \
	cicada_extract_scfg \
	cicada_extract_scfg_mpi \
	cicada_extract_score \
	cicada_extract_score_mpi \
	cicada_extract_tree \
	cicada_extract_tree_mpi \
	cicada_filter_cabocha \
	cicada_filter_extract \
	cicada_filter_extract_ghkm \
	cicada_filter_extract_phrase \
	cicada_filter_extract_scfg \
	cicada_filter_penntreebank \
	cicada_index_cluster \
	cicada_index_grammar \
	cicada_index_global_lexicon \
	cicada_index_tree_grammar \
	cicada_learn \
	cicada_learn_mpi \
	cicada_learn_online \
	cicada_learn_online_mpi \
	cicada_mert \
	cicada_mert_mpi \
	cicada_maxlike \
	cicada_maxlike_mpi \
	cicada_oracle \
	cicada_oracle_mpi \
	cicada_unite_forest \
	cicada_unite_lattice \
	cicada_unite_sentence \
	\
	mpish \
	mpipe

cicada_SOURCES = cicada.cpp
cicada_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mpi_SOURCES  = cicada_mpi.cpp
cicada_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_eval_SOURCES = cicada_eval.cpp
cicada_eval_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_ghkm_SOURCES = cicada_extract_ghkm.cpp
cicada_extract_ghkm_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_ghkm_mpi_SOURCES  = cicada_extract_ghkm_mpi.cpp
cicada_extract_ghkm_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_ghkm_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_phrase_SOURCES = cicada_extract_phrase.cpp
cicada_extract_phrase_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_phrase_mpi_SOURCES  = cicada_extract_phrase_mpi.cpp
cicada_extract_phrase_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_phrase_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_scfg_SOURCES = cicada_extract_scfg.cpp
cicada_extract_scfg_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_scfg_mpi_SOURCES  = cicada_extract_scfg_mpi.cpp
cicada_extract_scfg_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_scfg_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_score_SOURCES = cicada_extract_score.cpp
cicada_extract_score_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_score_mpi_SOURCES  = cicada_extract_score_mpi.cpp
cicada_extract_score_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_score_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_tree_SOURCES = cicada_extract_tree.cpp
cicada_extract_tree_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_tree_mpi_SOURCES  = cicada_extract_tree_mpi.cpp
cicada_extract_tree_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_tree_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_penntreebank_SOURCES = cicada_filter_penntreebank.cpp
cicada_filter_penntreebank_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_cabocha_SOURCES = cicada_filter_cabocha.cpp
cicada_filter_cabocha_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_extract_SOURCES = cicada_filter_extract.cpp
cicada_filter_extract_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_extract_ghkm_SOURCES = cicada_filter_extract_ghkm.cpp
cicada_filter_extract_ghkm_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_extract_phrase_SOURCES = cicada_filter_extract_phrase.cpp
cicada_filter_extract_phrase_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_extract_scfg_SOURCES = cicada_filter_extract_scfg.cpp
cicada_filter_extract_scfg_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_index_cluster_SOURCES = cicada_index_cluster.cpp
cicada_index_cluster_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_index_grammar_SOURCES = cicada_index_grammar.cpp
cicada_index_grammar_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_index_global_lexicon_SOURCES = cicada_index_global_lexicon.cpp
cicada_index_global_lexicon_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_index_tree_grammar_SOURCES = cicada_index_tree_grammar.cpp
cicada_index_tree_grammar_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_SOURCES = cicada_learn.cpp lbfgs.c
cicada_learn_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_mpi_SOURCES  = cicada_learn_mpi.cpp lbfgs.c
cicada_learn_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_learn_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_online_SOURCES = cicada_learn_online.cpp
cicada_learn_online_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_online_mpi_SOURCES  = cicada_learn_online_mpi.cpp
cicada_learn_online_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_learn_online_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mert_SOURCES = cicada_mert.cpp
cicada_mert_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mert_mpi_SOURCES  = cicada_mert_mpi.cpp
cicada_mert_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_mert_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_maxlike_SOURCES = cicada_maxlike.cpp lbfgs.c
cicada_maxlike_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_maxlike_mpi_SOURCES  = cicada_maxlike_mpi.cpp lbfgs.c
cicada_maxlike_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_maxlike_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_unite_forest_SOURCES = cicada_unite_forest.cpp
cicada_unite_forest_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_unite_lattice_SOURCES = cicada_unite_lattice.cpp
cicada_unite_lattice_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_unite_sentence_SOURCES = cicada_unite_sentence.cpp
cicada_unite_sentence_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_oracle_SOURCES = cicada_oracle.cpp
cicada_oracle_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_oracle_mpi_SOURCES  = cicada_oracle_mpi.cpp
cicada_oracle_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_oracle_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

mpish_SOURCES = mpish.cpp
mpish_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
mpish_LDADD   = $(MPI_LDFLAGS) $(boost_LDADD) $(perftools_LDADD)

mpipe_SOURCES = mpipe.cpp
mpipe_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
mpipe_LDADD   = $(MPI_LDFLAGS) $(boost_LDADD) $(perftools_LDADD)
