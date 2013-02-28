
if WITH_RST2MAN

man1_MANS = \
cicada.man \
cicada_mpi.man \
cicada_filter_alignment.man \
cicada_filter_bitext.man \
cicada_filter_chardet.man \
cicada_filter_charniak.man \
cicada_filter_config.man \
cicada_filter_config_moses.man \
cicada_filter_dependency.man

.rst.man: $*.rst
	$(RST2MAN) < $*.rst > $*.man

endif WITH_RST2MAN
