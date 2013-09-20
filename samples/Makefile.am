generated_models = \
scfg/grammar.bin \
s2t/grammar.bin \
t2s/grammar.bin

scfg/grammar.bin: scfg/grammar.bz2 ../progs/cicada_index_grammar
	../progs/cicada_index_grammar --input scfg/grammar.bz2 --output scfg/grammar.bin

s2t/grammar.bin: s2t/grammar.bz2 ../progs/cicada_index_tree_grammar
	../progs/cicada_index_tree_grammar --input s2t/grammar.bz2:cky=true --output s2t/grammar.bin

t2s/grammar.bin: t2s/grammar.bz2 ../progs/cicada_index_tree_grammar
	../progs/cicada_index_tree_grammar --input t2s/grammar.bz2 --output t2s/grammar.bin

noinst_DATA = $(generated_models)

sample_scfg_models = \
scfg/grammar.bz2 \
scfg/input.txt \
scfg/ngram.bin \
scfg/ngram.bz2 \
scfg/sample.sh \
scfg/weights


sample_s2t_models = \
s2t/grammar.bz2 \
s2t/sample.sh \
s2t/weights

sample_t2s_models = \
t2s/grammar.bz2 \
t2s/input.txt \
t2s/sample.sh \
t2s/weights

sample_models = \
$(sample_scfg_models) \
$(sample_s2t_models) \
$(sample_t2s_models)

sample_trains = \
kftt.30k

EXTRA_DIST = \
$(sample_models) \
$(sample_trains) \
cicada.config.hiero \
cicada.config.hiero-mbr \
cicada.config.string-to-tree \
cicada.config.tree-to-string \
dates_en.txt \
dates_zh_CN.txt \
numbers.en.txt \
numbers.ja.txt \
numbers.zh.txt \
rule.zone-wall \
weights.mbr
