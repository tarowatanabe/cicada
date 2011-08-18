AM_CPPFLAGS = -I.. -Iutils $(BOOST_CPPFLAGS)

noinst_LTLIBRARIES = libwn.la
noinst_HEADERS = wn.h wordnet.hpp

libwn_la_SOURCES = \
binsrch.c  \
morph.c    \
search.c   \
wnglobal.c \
wnhelp.c   \
wnrtl.c    \
wnutil.c \
wordnet.cpp

libwn_la_LIBADD = \
	$(BOOST_THREAD_LDFLAGS) $(BOOST_THREAD_LIBS) \
	$(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_SYSTEM_LDFLAGS) $(BOOST_SYSTEM_LIBS)

noinst_PROGRAMS = \
	wn_main

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done

wn_main_SOURCES = wn_main.cpp
wn_main_LDADD   = libwn.la

EXTRA_DIST = license.txt