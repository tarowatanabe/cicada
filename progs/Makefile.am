
AM_CPPFLAGS = -I.. $(ICU_CPPFLAGS) $(BOOST_CPPFLAGS) $(PTHREAD_FLAGS) $(TCMALLOC_CPPFLAGS) $(PROFILER_CPPFLAGS) $(SPARSEHASH_CPPFLAGS)

perftools_LDADD = $(JEMALLOC_LDFLAGS) $(TCMALLOC_LDFLAGS) $(PROFILER_LDFLAGS)

boost_LDADD = \
	$(BOOST_THREAD_LDFLAGS) $(BOOST_THREAD_LIBS) \
	$(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
	$(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_IOSTREAMS_LIBS) \
	$(BOOST_SYSTEM_LDFLAGS) $(BOOST_SYSTEM_LIBS)

LIBCICADA = $(top_builddir)/cicada/libcicada.la
LIBUTILS  = $(top_builddir)/utils/libutils.la
LIBLINEAR = $(top_builddir)/liblinear/liblinear.la

noinst_PROGRAMS = \
cicada_kbest_main \
cicada_text_main

cicada_kbest_main_SOURCES = cicada_kbest_main.cpp cicada_kbest_impl.hpp
cicada_kbest_main_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_text_main_SOURCES = cicada_text_main.cpp cicada_text_impl.hpp
cicada_text_main_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

bin_PROGRAMS = \
	cicada \
	cicada_mpi \
	cicada_alignment \
	cicada_cluster_word \
	cicada_convex_hull \
	cicada_convex_hull_mpi \
	cicada_convex_hull_kbest \
	cicada_convex_hull_kbest_mpi \
	cicada_eval \
	cicada_extract_ghkm \
	cicada_extract_ghkm_mpi \
	cicada_extract_phrase \
	cicada_extract_phrase_mpi \
	cicada_extract_rule \
	cicada_extract_scfg \
	cicada_extract_scfg_mpi \
	cicada_extract_score \
	cicada_extract_score_mpi \
	cicada_extract_sort \
	cicada_extract_tree \
	cicada_extract_tree_mpi \
	cicada_filter_bitext \
	cicada_filter_cabocha \
	cicada_filter_config \
	cicada_filter_extract \
	cicada_filter_extract_ghkm \
	cicada_filter_extract_phrase \
	cicada_filter_extract_scfg \
	cicada_filter_kbest \
	cicada_filter_charniak \
	cicada_filter_conll \
	cicada_filter_dependency \
	cicada_filter_penntreebank \
	cicada_filter_pialign \
	cicada_filter_terminal \
	cicada_filter_text \
	cicada_filter_weights \
	cicada_format_rbnf \
	cicada_grammar_coarse \
	cicada_grammar_learn \
	cicada_index_cluster \
	cicada_index_grammar \
	cicada_index_global_lexicon \
	cicada_index_tree_grammar \
	cicada_learn \
	cicada_learn_mpi \
	cicada_learn_kbest \
	cicada_learn_kbest_mpi \
	cicada_learn_online \
	cicada_learn_online_mpi \
	cicada_lexicon \
	cicada_lexicon_dice \
	cicada_lexicon_global \
	cicada_lexicon_global_mpi \
	cicada_lexicon_hmm \
	cicada_lexicon_model1 \
	cicada_maxlike \
	cicada_maxlike_mpi \
	cicada_mert \
	cicada_mert_mpi \
	cicada_mert_kbest \
	cicada_mert_kbest_mpi \
	cicada_oracle \
	cicada_oracle_mpi \
	cicada_oracle_kbest \
	cicada_oracle_kbest_mpi \
	cicada_unite_forest \
	cicada_unite_lattice \
	cicada_unite_sentence \
	\
	mpimap \
	mpish \
	mpipe \
	thrsh

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done

LBFGS_SOURCE = lbfgs.c arithmetic_ansi.h  arithmetic_sse_double.h  arithmetic_sse_float.h  lbfgs.h

cicada_SOURCES = cicada.cpp cicada_impl.hpp
cicada_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mpi_SOURCES  = cicada_mpi.cpp cicada_impl.hpp
cicada_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_alignment_SOURCES = cicada_alignment.cpp itg_alignment.hpp kuhn_munkres.hpp
cicada_alignment_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_cluster_word_SOURCES = cicada_cluster_word.cpp
cicada_cluster_word_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_convex_hull_SOURCES = cicada_convex_hull.cpp
cicada_convex_hull_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_convex_hull_mpi_SOURCES = cicada_convex_hull_mpi.cpp
cicada_convex_hull_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_convex_hull_mpi_LDADD   = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_convex_hull_kbest_SOURCES = cicada_convex_hull_kbest.cpp
cicada_convex_hull_kbest_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_convex_hull_kbest_mpi_SOURCES = cicada_convex_hull_kbest_mpi.cpp
cicada_convex_hull_kbest_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_convex_hull_kbest_mpi_LDADD   = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_eval_SOURCES = cicada_eval.cpp cicada_text_impl.hpp
cicada_eval_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_ghkm_SOURCES = cicada_extract_ghkm.cpp cicada_extract_ghkm_impl.hpp
cicada_extract_ghkm_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_ghkm_mpi_SOURCES  = cicada_extract_ghkm_mpi.cpp cicada_extract_ghkm_impl.hpp
cicada_extract_ghkm_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_ghkm_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_phrase_SOURCES = cicada_extract_phrase.cpp cicada_extract_phrase_impl.hpp
cicada_extract_phrase_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_phrase_mpi_SOURCES  = cicada_extract_phrase_mpi.cpp cicada_extract_phrase_impl.hpp
cicada_extract_phrase_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_phrase_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_rule_SOURCES = cicada_extract_rule.cpp
cicada_extract_rule_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_scfg_SOURCES = cicada_extract_scfg.cpp cicada_extract_scfg_impl.hpp
cicada_extract_scfg_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_scfg_mpi_SOURCES  = cicada_extract_scfg_mpi.cpp cicada_extract_scfg_impl.hpp
cicada_extract_scfg_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_scfg_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_score_SOURCES = cicada_extract_score.cpp cicada_extract_score_impl.hpp
cicada_extract_score_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_score_mpi_SOURCES  = cicada_extract_score_mpi.cpp cicada_extract_score_impl.hpp
cicada_extract_score_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_score_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_sort_SOURCES = cicada_extract_sort.cpp cicada_extract_score_impl.hpp
cicada_extract_sort_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_tree_SOURCES = cicada_extract_tree.cpp cicada_extract_tree_impl.hpp
cicada_extract_tree_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_extract_tree_mpi_SOURCES  = cicada_extract_tree_mpi.cpp cicada_extract_tree_impl.hpp
cicada_extract_tree_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_extract_tree_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_charniak_SOURCES = cicada_filter_charniak.cpp
cicada_filter_charniak_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_conll_SOURCES = cicada_filter_conll.cpp
cicada_filter_conll_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_dependency_SOURCES = cicada_filter_dependency.cpp
cicada_filter_dependency_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_penntreebank_SOURCES = cicada_filter_penntreebank.cpp
cicada_filter_penntreebank_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_pialign_SOURCES = cicada_filter_pialign.cpp
cicada_filter_pialign_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_bitext_SOURCES = cicada_filter_bitext.cpp
cicada_filter_bitext_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_text_SOURCES = cicada_filter_text.cpp
cicada_filter_text_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_cabocha_SOURCES = cicada_filter_cabocha.cpp
cicada_filter_cabocha_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_config_SOURCES = cicada_filter_config.cpp
cicada_filter_config_LDADD   = $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_terminal_SOURCES = cicada_filter_terminal.cpp
cicada_filter_terminal_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_extract_SOURCES = cicada_filter_extract.cpp cicada_filter_extract_impl.hpp
cicada_filter_extract_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_extract_ghkm_SOURCES = cicada_filter_extract_ghkm.cpp cicada_filter_extract_impl.hpp
cicada_filter_extract_ghkm_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_extract_phrase_SOURCES = cicada_filter_extract_phrase.cpp cicada_filter_extract_impl.hpp
cicada_filter_extract_phrase_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_extract_scfg_SOURCES = cicada_filter_extract_scfg.cpp cicada_filter_extract_impl.hpp
cicada_filter_extract_scfg_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_kbest_SOURCES = cicada_filter_kbest.cpp
cicada_filter_kbest_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_filter_weights_SOURCES = cicada_filter_weights.cpp
cicada_filter_weights_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_format_rbnf_SOURCES = cicada_format_rbnf.cpp
cicada_format_rbnf_LDADD   = $(LIBUTILS) $(ICU_LDFLAGS) $(boost_LDADD) $(perftools_LDADD)

cicada_grammar_coarse_SOURCES = cicada_grammar_coarse.cpp
cicada_grammar_coarse_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_grammar_learn_SOURCES  = cicada_grammar_learn.cpp
cicada_grammar_learn_LDADD    = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_index_cluster_SOURCES = cicada_index_cluster.cpp
cicada_index_cluster_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_index_grammar_SOURCES = cicada_index_grammar.cpp
cicada_index_grammar_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_index_global_lexicon_SOURCES = cicada_index_global_lexicon.cpp
cicada_index_global_lexicon_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_index_tree_grammar_SOURCES = cicada_index_tree_grammar.cpp
cicada_index_tree_grammar_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_SOURCES = cicada_learn.cpp $(LBFGS_SOURCE) cicada_learn_impl.hpp
cicada_learn_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_mpi_SOURCES  = cicada_learn_mpi.cpp $(LBFGS_SOURCE) cicada_learn_impl.hpp
cicada_learn_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_learn_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_kbest_SOURCES = cicada_learn_kbest.cpp $(LBFGS_SOURCE) cicada_kbest_impl.hpp
cicada_learn_kbest_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD) $(LIBLINEAR)

cicada_learn_kbest_mpi_SOURCES  = cicada_learn_kbest_mpi.cpp $(LBFGS_SOURCE) cicada_kbest_impl.hpp
cicada_learn_kbest_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_learn_kbest_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_online_SOURCES = cicada_learn_online.cpp cicada_learn_online_impl.hpp
cicada_learn_online_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_learn_online_mpi_SOURCES  = cicada_learn_online_mpi.cpp cicada_learn_online_impl.hpp
cicada_learn_online_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_learn_online_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_lexicon_SOURCES = cicada_lexicon.cpp cicada_lexicon_impl.hpp
cicada_lexicon_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_lexicon_dice_SOURCES = cicada_lexicon_dice.cpp
cicada_lexicon_dice_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_lexicon_global_SOURCES = cicada_lexicon_global.cpp $(LBFGS_SOURCE) cicada_lexicon_global_impl.hpp
cicada_lexicon_global_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_lexicon_global_mpi_SOURCES = cicada_lexicon_global_mpi.cpp $(LBFGS_SOURCE) cicada_lexicon_global_impl.hpp
cicada_lexicon_global_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_lexicon_global_mpi_LDADD   = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_lexicon_hmm_SOURCES = cicada_lexicon_hmm.cpp cicada_lexicon_hmm_impl.hpp cicada_lexicon_model1_impl.hpp cicada_lexicon_maximize_impl.hpp kuhn_munkres.hpp
cicada_lexicon_hmm_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_lexicon_model1_SOURCES = cicada_lexicon_model1.cpp cicada_lexicon_model1_impl.hpp cicada_lexicon_maximize_impl.hpp kuhn_munkres.hpp
cicada_lexicon_model1_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mert_SOURCES = cicada_mert.cpp cicada_text_impl.hpp
cicada_mert_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mert_mpi_SOURCES  = cicada_mert_mpi.cpp cicada_text_impl.hpp
cicada_mert_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_mert_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mert_kbest_SOURCES = cicada_mert_kbest.cpp cicada_text_impl.hpp
cicada_mert_kbest_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_mert_kbest_mpi_SOURCES  = cicada_mert_kbest_mpi.cpp cicada_text_impl.hpp
cicada_mert_kbest_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_mert_kbest_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_maxlike_SOURCES = cicada_maxlike.cpp $(LBFGS_SOURCE) cicada_maxlike_impl.hpp cicada_text_impl.hpp
cicada_maxlike_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_maxlike_mpi_SOURCES  = cicada_maxlike_mpi.cpp $(LBFGS_SOURCE) cicada_maxlike_impl.hpp cicada_text_impl.hpp
cicada_maxlike_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_maxlike_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_unite_forest_SOURCES = cicada_unite_forest.cpp
cicada_unite_forest_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_unite_lattice_SOURCES = cicada_unite_lattice.cpp
cicada_unite_lattice_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_unite_sentence_SOURCES = cicada_unite_sentence.cpp
cicada_unite_sentence_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_oracle_SOURCES = cicada_oracle.cpp cicada_text_impl.hpp
cicada_oracle_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_oracle_mpi_SOURCES  = cicada_oracle_mpi.cpp cicada_text_impl.hpp
cicada_oracle_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_oracle_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_oracle_kbest_SOURCES = cicada_oracle_kbest.cpp cicada_text_impl.hpp
cicada_oracle_kbest_LDADD   = $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

cicada_oracle_kbest_mpi_SOURCES  = cicada_oracle_kbest_mpi.cpp cicada_text_impl.hpp
cicada_oracle_kbest_mpi_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
cicada_oracle_kbest_mpi_LDADD    = $(MPI_LDFLAGS) $(LIBCICADA) $(LIBUTILS) $(boost_LDADD) $(perftools_LDADD)

mpimap_SOURCES = mpimap.cpp
mpimap_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
mpimap_LDADD   = $(MPI_LDFLAGS) $(boost_LDADD) $(perftools_LDADD)

mpish_SOURCES = mpish.cpp
mpish_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
mpish_LDADD   = $(MPI_LDFLAGS) $(boost_LDADD) $(perftools_LDADD)

mpipe_SOURCES = mpipe.cpp
mpipe_CPPFLAGS = $(MPI_CPPFLAGS) $(AM_CPPFLAGS)
mpipe_LDADD   = $(MPI_LDFLAGS) $(boost_LDADD) $(perftools_LDADD)

thrsh_SOURCES = thrsh.cpp
thrsh_CPPFLAGS = $(AM_CPPFLAGS)
thrsh_LDADD   = $(boost_LDADD) $(perftools_LDADD)
