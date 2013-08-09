generated_models = \
scfg/grammar.bin \
s2t/grammar.bin \
t2s/grammar.bin

scfg/grammar.bin: scfg/grammar.bz2
	../progs/cicada_index_grammar --input scfg/grammar.bz2 --output scfg/grammar.bin

s2t/grammar.bin: s2t/grammar.bz2
	../progs/cicada_index_tree_grammar --input s2t/grammar.bz2:cky=true --output s2t/grammar.bin

t2s/grammar.bin: t2s/grammar.bz2
	../progs/cicada_index_tree_grammar --input t2s/grammar.bz2 --output t2s/grammar.bin

noinst_DATA = $(generated_models)

sample_models = \
scfg \
s2t \
t2s


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
