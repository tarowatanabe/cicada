
if WITH_RST2MAN

man1_MANS = \
cicada.man \
cicada_mpi.man

.rst.man: $*.rst
	$(RST2MAN) < $*.rst > $*.man

endif WITH_RST2MAN
