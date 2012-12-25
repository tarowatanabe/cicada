AM_CPPFLAGS = -I.. -I../eigen $(BOOST_CPPFLAGS) $(ICU_CPPFLAGS) $(PTHREAD_FLAGS) $(MSGPACK_CPPFLAGS)

LIBUTILS=$(top_builddir)/utils/libutils.la

lib_LTLIBRARIES = libcicada.la

noinst_HEADERS =\
msgpack_main_impl.hpp

noinst_LTLIBRARIES = \
libcicada-feature.la \
libcicada-matcher.la \
libcicada-neuron.la \
libcicada-semiring.la \
libcicada-eval.la \
libcicada-format.la \
libcicada-head.la \
libcicada-optimize.la \
libcicada-operation.la \
libcicada-signature.la  \
libcicada-stemmer.la 

libcicadaincludedir= $(pkgincludedir)/cicada

evaldir      = $(pkgincludedir)/cicada/eval
featuredir   = $(pkgincludedir)/cicada/feature
formatdir    = $(pkgincludedir)/cicada/format
headdir      = $(pkgincludedir)/cicada/head
matcherdir   = $(pkgincludedir)/cicada/matcher
msgpackdir   = $(pkgincludedir)/cicada/msgpack
neurondir    = $(pkgincludedir)/cicada/neuron
operationdir = $(pkgincludedir)/cicada/operation
optimizedir  = $(pkgincludedir)/cicada/optimize
semiringdir  = $(pkgincludedir)/cicada/semiring
signaturedir = $(pkgincludedir)/cicada/signature
stemmerdir   = $(pkgincludedir)/cicada/stemmer
tokenizerdir = $(pkgincludedir)/cicada/tokenizer

dist_semiring_HEADERS = \
semiring/envelope.hpp \
semiring/expectation.hpp \
semiring/log.hpp \
semiring/logprob.hpp \
semiring/pair.hpp \
semiring/tropical.hpp \
semiring/traits.hpp \
semiring/tuple.hpp

libcicada_semiring_la_SOURCES = \
semiring/envelope.cpp
libcicada_semiring_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_eval_HEADERS = \
eval/combined.hpp \
eval/f.hpp \
eval/bleu.hpp \
eval/bleus.hpp \
eval/cder.hpp \
eval/depeval.hpp \
eval/inv_wer.hpp \
eval/parseval.hpp \
eval/per.hpp \
eval/ribes.hpp \
eval/score.hpp \
eval/sb.hpp \
eval/sk.hpp \
eval/ter.hpp \
eval/wer.hpp \
eval/wlcs.hpp

libcicada_eval_la_SOURCES = \
eval/bleu.cpp \
eval/bleus.cpp \
eval/cder.cpp \
eval/combined.cpp \
eval/decode.hpp \
eval/depeval.cpp \
eval/encode.hpp \
eval/f.cpp \
eval/inv_wer.cpp \
eval/parseval.cpp \
eval/per.cpp \
eval/ribes.cpp \
eval/score.cpp \
eval/sb.cpp \
eval/sk.cpp \
eval/ter.cpp \
eval/wer.cpp \
eval/wlcs.cpp
libcicada_eval_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_format_HEADERS = \
format/date.hpp \
format/number.hpp

libcicada_format_la_SOURCES = \
format/date.cpp \
format/number.cpp

libcicada_format_la_LIBADD = $(ICU_LDFLAGS)
libcicada_format_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_head_HEADERS = \
head/chinese.hpp \
head/collins.hpp \
head/collins_modified.hpp 

libcicada_head_la_SOURCES = \
head/chinese.cpp \
head/collins.cpp \
head/collins_modified.cpp

libcicada_head_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_optimize_HEADERS = \
optimize/edmonds_optimum_branching.hpp \
optimize/edmonds_optimum_branching_impl.hpp \
optimize/line_search.hpp \
optimize/mst.hpp \
optimize/powell.hpp \
optimize/qp_dcd.hpp \
optimize/qp_simplex.hpp \
optimize/qp_smo.hpp

libcicada_optimize_la_SOURCES = \
optimize/line_search.cpp
libcicada_optimize_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_feature_HEADERS = \
feature/alignment.hpp \
feature/antecedent.hpp \
feature/bleu.hpp \
feature/bleus.hpp \
feature/bleu_expected.hpp \
feature/bleu_linear.hpp \
feature/bleu_multi.hpp \
feature/dependency.hpp \
feature/depeval.hpp \
feature/deletion.hpp \
feature/distortion.hpp \
feature/feature_builder.hpp \
feature/global_lexicon.hpp \
feature/insertion.hpp \
feature/lexicalized_reordering.hpp \
feature/lexicon.hpp \
feature/neighbours.hpp \
feature/ngram.hpp \
feature/ngram_pyp.hpp \
feature/ngram_tree.hpp \
feature/parent.hpp \
feature/penalty.hpp \
feature/permute.hpp \
feature/rule_shape.hpp \
feature/scorer.hpp \
feature/sgml_tag.hpp \
feature/span.hpp \
feature/sparse_lexicon.hpp \
feature/sparse_ngram.hpp \
feature/variational.hpp \
feature/vocabulary.hpp

libcicada_feature_la_SOURCES = \
feature/alignment.cpp \
feature/antecedent.cpp \
feature/bleu.cpp \
feature/bleus.cpp \
feature/bleu_expected.cpp \
feature/bleu_linear.cpp \
feature/bleu_multi.cpp \
feature/distortion.cpp \
feature/dependency.cpp \
feature/depeval.cpp \
feature/deletion.cpp \
feature/global_lexicon.cpp \
feature/insertion.cpp \
feature/lexicalized_reordering.cpp \
feature/lexicon.cpp \
feature/neighbours.cpp \
feature/ngram.cpp \
feature/ngram_pyp.cpp \
feature/ngram_tree.cpp \
feature/parent.cpp \
feature/penalty.cpp \
feature/permute.cpp \
feature/rule_shape.cpp \
feature/sgml_tag.cpp \
feature/span.cpp \
feature/sparse_lexicon.cpp \
feature/sparse_ngram.cpp \
feature/variational.cpp \
feature/vocabulary.cpp

libcicada_feature_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_signature_HEADERS = \
signature/chinese.hpp  \
signature/english.hpp 

libcicada_signature_la_SOURCES = \
signature/chinese.cpp \
signature/english.cpp

libcicada_signature_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_stemmer_HEADERS = \
stemmer/arabic.hpp \
stemmer/prefix.hpp \
stemmer/suffix.hpp \
stemmer/katakana.hpp \
stemmer/latin.hpp \
stemmer/lower.hpp \
stemmer/upper.hpp \
stemmer/fullwidth.hpp \
stemmer/halfwidth.hpp \
stemmer/simplified.hpp \
stemmer/traditional.hpp \
stemmer/digit.hpp \
stemmer/nfkc.hpp \
stemmer/nfkd.hpp \
stemmer/snowball.hpp

libcicada_stemmer_la_SOURCES = \
stemmer/arabic.cpp \
stemmer/prefix.cpp \
stemmer/suffix.cpp \
stemmer/katakana.cpp \
stemmer/latin.cpp \
stemmer/lower.cpp \
stemmer/upper.cpp \
stemmer/fullwidth.cpp \
stemmer/halfwidth.cpp \
stemmer/simplified.cpp \
stemmer/traditional.cpp \
stemmer/digit.cpp \
stemmer/nfkc.cpp \
stemmer/nfkd.cpp \
stemmer/snowball.cpp

libcicada_stemmer_la_CPPFLAGS = $(AM_CPPFLAGS)
libcicada_stemmer_la_LIBADD = $(top_builddir)/libstemmer_c/libstemmer.la

dist_tokenizer_HEADERS = \
tokenizer/characters.hpp \
tokenizer/lower.hpp \
tokenizer/nist.hpp \
tokenizer/nonascii.hpp \
tokenizer/penntreebank.hpp \
tokenizer/stemmer.hpp \
tokenizer/tokenize.hpp

dist_matcher_HEADERS = \
matcher/lower.hpp \
matcher/stemmer.hpp \
matcher/wordnet.hpp 

libcicada_matcher_la_SOURCES = \
matcher/wordnet.cpp

libcicada_matcher_la_CPPFLAGS = $(AM_CPPFLAGS)
libcicada_matcher_la_LIBADD = $(top_builddir)/wn/libwn.la

dist_msgpack_HEADERS = \
msgpack/alignment.hpp \
msgpack/attribute.hpp \
msgpack/attribute_vector.hpp \
msgpack/dependency.hpp \
msgpack/feature.hpp \
msgpack/feature_vector.hpp \
msgpack/feature_vector_linear.hpp \
msgpack/hypergraph.hpp \
msgpack/lattice.hpp \
msgpack/ngram_count_set.hpp \
msgpack/operation.hpp \
msgpack/rule.hpp \
msgpack/sentence.hpp \
msgpack/sentence_vector.hpp \
msgpack/span_vector.hpp \
msgpack/statistics.hpp \
msgpack/symbol.hpp \
msgpack/symbol_vector.hpp \
msgpack/tree_rule.hpp \
msgpack/weight_vector.hpp

dist_neuron_HEADERS = \
neuron/abs.hpp \
neuron/concat.hpp \
neuron/convolution.hpp \
neuron/copy.hpp \
neuron/exp.hpp \
neuron/hardtanh.hpp \
neuron/layer.hpp \
neuron/linear.hpp \
neuron/log.hpp \
neuron/lookup.hpp \
neuron/max.hpp \
neuron/mean.hpp \
neuron/min.hpp \
neuron/neuron.hpp \
neuron/parallel.hpp \
neuron/power.hpp \
neuron/sequential.hpp \
neuron/sigmoid.hpp \
neuron/softmax.hpp \
neuron/sqrt.hpp \
neuron/square.hpp \
neuron/sum.hpp \
neuron/tanh.hpp

libcicada_neuron_la_SOURCES = \
neuron/abs.cpp \
neuron/concat.cpp \
neuron/convolution.cpp \
neuron/copy.cpp \
neuron/hardtanh.cpp \
neuron/exp.cpp \
neuron/linear.cpp \
neuron/log.cpp \
neuron/lookup.cpp \
neuron/max.cpp \
neuron/mean.cpp \
neuron/min.cpp \
neuron/parallel.cpp \
neuron/power.cpp \
neuron/sequential.cpp \
neuron/sigmoid.cpp \
neuron/softmax.cpp \
neuron/sqrt.cpp \
neuron/square.cpp \
neuron/sum.cpp \
neuron/tanh.cpp

libcicada_neuron_la_CPPFLAGS = $(AM_CPPFLAGS)

dist_operation_HEADERS = \
operation/apply.hpp \
operation/binarize.hpp \
operation/clear.hpp \
operation/compose.hpp \
operation/debinarize.hpp \
operation/expand_ngram.hpp \
operation/expected_ngram.hpp \
operation/functional.hpp \
operation/generate.hpp \
operation/intersect.hpp \
operation/normalize.hpp \
operation/output.hpp \
operation/parse.hpp \
operation/permute.hpp \
operation/posterior.hpp \
operation/prune.hpp \
operation/push_bos_eos.hpp \
operation/push_head.hpp \
operation/push_weights.hpp \
operation/remove_annotation.hpp \
operation/remove_bos_eos.hpp \
operation/remove_epsilon.hpp \
operation/remove_feature.hpp \
operation/remove_head.hpp \
operation/remove_sgml_tag.hpp \
operation/remove_unary.hpp \
operation/sort_tail.hpp \
operation/sort_topologically.hpp \
operation/span_forest.hpp \
operation/traversal.hpp \
operation/viterbi.hpp

libcicada_operation_la_SOURCES = \
operation/apply.cpp \
operation/binarize.cpp \
operation/clear.cpp \
operation/compose.cpp \
operation/debinarize.cpp \
operation/expand_ngram.cpp \
operation/expected_ngram.cpp \
operation/intersect.cpp \
operation/generate.cpp \
operation/normalize.cpp \
operation/output.cpp \
operation/parse.cpp \
operation/permute.cpp \
operation/posterior.cpp \
operation/prune.cpp \
operation/push_bos_eos.cpp \
operation/push_head.cpp \
operation/push_weights.cpp \
operation/remove_annotation.cpp \
operation/remove_bos_eos.cpp \
operation/remove_epsilon.cpp \
operation/remove_feature.cpp \
operation/remove_head.cpp \
operation/remove_sgml_tag.cpp \
operation/remove_unary.cpp \
operation/sort_tail.cpp \
operation/sort_topologically.cpp \
operation/span_forest.cpp \
operation/viterbi.cpp
libcicada_operation_la_CPPFLAGS = $(AM_CPPFLAGS)


libcicadainclude_HEADERS = \
alignment.hpp \
apply.hpp \
apply_cube_grow.hpp \
apply_cube_grow_coarse.hpp \
apply_cube_prune.hpp \
apply_exact.hpp \
apply_incremental.hpp \
apply_state_less.hpp \
attribute.hpp \
attribute_vector.hpp \
attribute_vector_codec.hpp \
binarize.hpp \
binarize_all.hpp \
binarize_base.hpp \
binarize_cyk.hpp \
binarize_dependency.hpp \
binarize_left.hpp \
binarize_right.hpp \
cluster.hpp \
cluster_stemmer.hpp \
compose.hpp \
compose_alignment.hpp \
compose_cky.hpp \
compose_dependency.hpp \
compose_dependency_arc_standard.hpp \
compose_dependency_arc_eager.hpp \
compose_dependency_degree2.hpp \
compose_dependency_hybrid.hpp \
compose_dependency_top_down.hpp \
compose_dependency_tree.hpp \
compose_earley.hpp \
compose_grammar.hpp \
compose_phrase.hpp \
compose_tree.hpp \
compose_tree_cky.hpp \
debinarize.hpp \
dependency.hpp \
discounter.hpp \
dot_product.hpp \
edit_distance_forest.hpp \
eval.hpp \
expected_ngram.hpp \
feature.hpp \
feature_function.hpp \
feature_vector.hpp \
feature_vector_codec.hpp \
feature_vector_compact.hpp \
feature_vector_linear.hpp \
format.hpp \
generate.hpp \
generate_earley.hpp \
global_lexicon.hpp \
grammar.hpp \
grammar_format.hpp \
grammar_mutable.hpp \
grammar_simple.hpp \
grammar_static.hpp \
grammar_shared.hpp \
grammar_unknown.hpp \
graphviz.hpp \
head_finder.hpp \
hypergraph.hpp \
inside_outside.hpp \
intersect.hpp \
kbest.hpp \
lattice.hpp \
lexicon.hpp \
matcher.hpp \
model.hpp \
ngram.hpp \
ngram_cache.hpp \
ngram_index.hpp \
ngram_count_set.hpp \
ngram_pyp.hpp \
operation.hpp \
operation_set.hpp \
optimize.hpp \
optimize_qp.hpp \
parameter.hpp \
parse.hpp \
parse_agenda.hpp \
parse_cky.hpp \
parse_coarse.hpp \
parse_phrase.hpp \
parse_tree.hpp \
parse_tree_cky.hpp \
permute.hpp \
permute_deterministic.hpp \
posterior.hpp \
prune.hpp \
prune_beam.hpp \
prune_density.hpp \
prune_edge.hpp \
prune_kbest.hpp \
push_bos_eos.hpp \
push_head.hpp \
push_weights.hpp \
push_weights_frontier.hpp \
push_weights_left.hpp \
push_weights_root.hpp \
quantizer.hpp \
query.hpp \
query_cky.hpp \
query_tree.hpp \
query_tree_cky.hpp \
remove_annotation.hpp \
remove_bos_eos.hpp \
remove_epsilon.hpp \
remove_head.hpp \
remove_non_terminal.hpp \
remove_sgml_tag.hpp \
remove_symbol.hpp \
remove_unary.hpp \
rule.hpp \
semiring.hpp \
sentence.hpp \
sentence_vector.hpp \
signature.hpp \
sort.hpp \
sort_tail.hpp \
sort_topologically.hpp \
span_node.hpp \
span_edge.hpp \
span_forest.hpp \
span_vector.hpp \
statistics.hpp \
stemmer.hpp \
tokenizer.hpp \
treebank.hpp \
symbol.hpp \
symbol_vector.hpp \
symbol_vector_compact.hpp \
transducer.hpp \
tree_grammar.hpp \
tree_grammar_mutable.hpp \
tree_grammar_static.hpp \
tree_grammar_shared.hpp \
tree_grammar_simple.hpp \
tree_rule.hpp \
tree_rule_codec.hpp \
tree_rule_compact.hpp \
tree_transducer.hpp \
unite.hpp \
viterbi.hpp \
vocab.hpp \
weight_vector.hpp

libcicada_la_SOURCES = \
alignment.cpp \
attribute.cpp \
attribute_vector.cpp \
cluster.cpp \
dependency.cpp \
feature.cpp \
feature_function.cpp \
format.cpp \
global_lexicon.cpp \
grammar_format.cpp \
grammar_mutable.cpp \
grammar_static.cpp \
grammar_unknown.cpp \
graphviz.cpp \
head_finder.cpp \
hypergraph.cpp \
lattice.cpp \
lexicon.cpp \
matcher.cpp \
model.cpp \
ngram.cpp \
ngram_index.cpp \
ngram_pyp.cpp \
operation.cpp \
operation_set.cpp \
parameter.cpp \
rule.cpp \
sentence.cpp \
sentence_vector.cpp \
signature.cpp \
span_vector.cpp \
statistics.cpp \
stemmer.cpp \
tokenizer.cpp \
treebank.cpp \
symbol.cpp \
symbol_vector.cpp \
transducer.cpp \
tree_grammar_mutable.cpp \
tree_grammar_static.cpp \
tree_rule.cpp \
tree_rule_compact.cpp \
tree_transducer.cpp \
vocab.cpp

libcicada_la_LDFLAGS = \
	$(PTHREAD_FLAGS) \
	$(BOOST_THREAD_LDFLAGS) \
	$(BOOST_FILESYSTEM_LDFLAGS) \
	$(BOOST_IOSTREAMS_LDFLAGS) \
	$(BOOST_REGEX_LDFLAGS) \
	-version-info $(CICADA_LTVERSION)

libcicada_la_LIBADD = \
	libcicada-eval.la \
	libcicada-format.la \
	libcicada-head.la \
	libcicada-feature.la \
	libcicada-matcher.la \
	libcicada-neuron.la \
	libcicada-optimize.la \
	libcicada-operation.la \
	libcicada-semiring.la \
	libcicada-signature.la \
	libcicada-stemmer.la \
	\
	$(LIBUTILS) \
	$(BOOST_THREAD_LIBS) \
	$(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_IOSTREAMS_LIBS) \
	$(BOOST_REGEX_LIBS) \
	$(ICU_LDFLAGS)

noinst_PROGRAMS = \
alignment_main \
attribute_vector_main \
cluster_main \
eval_main \
feature_vector_main \
format_main \
grammar_mutable_main \
grammar_static_main \
grammar_unknown_main \
hypergraph_main \
lattice_main \
lexicon_main \
matcher_main \
neuron_main \
ngram_pyp_main \
optimize_qp_main \
parameter_main \
rule_main \
sentence_main \
sentence_vector_main \
signature_main \
span_vector_main \
stemmer_main \
symbol_main \
symbol_vector_main \
tokenizer_main \
transducer_main \
tree_grammar_mutable_main \
tree_grammar_static_main \
tree_rule_main

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done

alignment_main_SOURCES = alignment_main.cpp
alignment_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

attribute_vector_main_SOURCES = attribute_vector_main.cpp
attribute_vector_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

cluster_main_SOURCES = cluster_main.cpp
cluster_main_LDADD = libcicada.la

eval_main_SOURCES = eval_main.cpp
eval_main_LDADD = libcicada.la

feature_vector_main_SOURCES = feature_vector_main.cpp symbol.cpp vocab.cpp feature.cpp
feature_vector_main_CPPFLAGS = $(AM_CPPFLAGS)
feature_vector_main_LDADD = $(LIBUTILS) $(MSGPACK_LDFLAGS)

format_main_SOURCES = format_main.cpp
format_main_LDADD = libcicada.la

grammar_mutable_main_SOURCES = grammar_mutable_main.cpp
grammar_mutable_main_LDADD = libcicada.la

grammar_static_main_SOURCES = grammar_static_main.cpp
grammar_static_main_LDADD = libcicada.la

grammar_unknown_main_SOURCES = grammar_unknown_main.cpp
grammar_unknown_main_LDADD = libcicada.la

hypergraph_main_SOURCES = hypergraph_main.cpp
hypergraph_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

lattice_main_SOURCES = lattice_main.cpp
lattice_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

lexicon_main_SOURCES = lexicon_main.cpp
lexicon_main_LDADD = libcicada.la

matcher_main_SOURCES = matcher_main.cpp
matcher_main_LDADD = libcicada.la

neuron_main_SOURCES = neuron_main.cpp
neuron_main_LDADD = libcicada.la

ngram_pyp_main_SOURCES = ngram_pyp_main.cpp
ngram_pyp_main_LDADD = libcicada.la

optimize_qp_main_SOURCES = optimize_qp_main.cpp
optimize_qp_main_LDADD = libcicada.la

parameter_main_SOURCES = parameter_main.cpp
parameter_main_LDADD = libcicada.la

rule_main_SOURCES = rule_main.cpp
rule_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

sentence_main_SOURCES = sentence_main.cpp
sentence_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

sentence_vector_main_SOURCES = sentence_vector_main.cpp
sentence_vector_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

signature_main_SOURCES = signature_main.cpp
signature_main_LDADD = libcicada.la

span_vector_main_SOURCES = span_vector_main.cpp
span_vector_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

stemmer_main_SOURCES = stemmer_main.cpp
stemmer_main_LDADD = libcicada.la

symbol_main_SOURCES = symbol_main.cpp
symbol_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

symbol_vector_main_SOURCES = symbol_vector_main.cpp
symbol_vector_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)

tokenizer_main_SOURCES = tokenizer_main.cpp
tokenizer_main_LDADD = libcicada.la

transducer_main_SOURCES = transducer_main.cpp
transducer_main_LDADD = libcicada.la

tree_grammar_mutable_main_SOURCES = tree_grammar_mutable_main.cpp
tree_grammar_mutable_main_LDADD = libcicada.la

tree_grammar_static_main_SOURCES = tree_grammar_static_main.cpp
tree_grammar_static_main_LDADD = libcicada.la

tree_rule_main_SOURCES = tree_rule_main.cpp
tree_rule_main_LDADD = libcicada.la $(MSGPACK_LDFLAGS)
