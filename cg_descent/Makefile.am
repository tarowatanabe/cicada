AM_CPPFLAGS = -I$(top_builddir)

LIBCG_DESCENT=libcg_descent.la

noinst_LTLIBRARIES = libcg_descent.la
noinst_HEADERS = \
	cg.hpp \
	cg_user.h \
	cg_descent.h \
	cg_blas.h

libcg_descent_la_SOURCES = \
	cg_descent.c

libcg_descent_la_LDFLAGS = -version-info $(CICADA_LTVERSION)

noinst_PROGRAMS = driver1 driver2 driver3 driver4 driver5

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

driver1_SOURCES = driver1.c
driver1_LDADD = libcg_descent.la

driver2_SOURCES = driver2.c
driver2_LDADD = libcg_descent.la

driver3_SOURCES = driver3.c
driver3_LDADD = libcg_descent.la

driver4_SOURCES = driver4.c
driver4_LDADD = libcg_descent.la

driver5_SOURCES = driver5.c
driver5_LDADD = libcg_descent.la

