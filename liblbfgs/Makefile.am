AM_CPPFLAGS = -I$(top_builddir)

LIBLBFGS=liblbfgs.la

noinst_LTLIBRARIES = liblbfgs.la
noinst_HEADERS = \
	arithmetic_ansi.h \
	arithmetic_sse_double.h \
	arithmetic_sse_float.h \
	lbfgs.h \
	lbfgs_error.hpp

liblbfgs_la_SOURCES = \
	lbfgs.c

liblbfgs_la_LDFLAGS = -version-info $(CICADA_LTVERSION)

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	echo "*.model" >> .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done

