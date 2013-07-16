
if WITH_RST2MAN

man1_MANS = \
cicada-alignment.py.man \
cicada-extract.py.man \
cicada.man \
cicada_mpi.man \
cicada_filter_alignment.man \
cicada_filter_bitext.man \
cicada_filter_chardet.man \
cicada_filter_charniak.man \
cicada_filter_config.man \
cicada_filter_config_moses.man \
cicada_filter_dependency.man

EXTRA_DIST = \
alignment.rst \
attributes.rst \
cicada-alignment.py.rst \
cicada-extract.py.rst \
cicada.rst \
cicada_filter_alignment.rst \
cicada_filter_bitext.rst \
cicada_filter_chardet.rst \
cicada_filter_charniak.rst \
cicada_filter_config.rst \
cicada_filter_config_moses.rst \
cicada_filter_dependency.rst \
cicada_mpi.rst \
concepts.rst \
eval.rst \
extract.rst \
features.rst \
grammar.rst \
hypergraph.rst \
index.rst \
lattice.rst \
learn.rst \
linguistic.rst \
mpi.rst \
non-terminals.rst \
operation.rst \
parameter.rst \
pyp.rst \
tree-grammar.rst

.rst.man: $*.rst
	$(RST2MAN) < $*.rst > $*.man

endif WITH_RST2MAN
