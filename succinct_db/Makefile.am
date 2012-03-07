AM_CPPFLAGS = -I.. $(BOOST_CPPFLAGS)

LIBUTILS=$(top_builddir)/utils/libutils.la

perftools_LDADD = $(JEMALLOC_LDFLAGS) $(TCMALLOC_LDFLAGS) $(PROFILER_LDFLAGS)

boost_LDADD = \
	$(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
	$(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_IOSTREAMS_LIBS) \
	$(BOOST_THREAD_LDFLAGS) $(BOOST_THREAD_LIBS)

succinct_dbincludedir = $(pkgincludedir)/succinct_db
succinct_dbinclude_HEADERS = \
	succinct_trie.hpp \
	succinct_hash.hpp \
	succinct_trie_db.hpp \
	succinct_trie_database.hpp \
	succinct_hash_db.hpp

noinst_PROGRAMS = \
	succinct_trie_main \
	succinct_hash_main \
	succinct_trie_db_main \
	succinct_trie_database_main

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done

succinct_trie_main_SOURCES = succinct_trie_main.cpp
succinct_trie_main_LDADD = $(boost_LDADD) $(LIBUTILS) $(perftools_LDADD)

succinct_hash_main_SOURCES = succinct_hash_main.cpp
succinct_hash_main_LDADD = $(boost_LDADD) $(LIBUTILS) $(perftools_LDADD)

succinct_trie_db_main_SOURCES = succinct_trie_db_main.cpp
succinct_trie_db_main_LDADD = $(boost_LDADD) $(LIBUTILS) $(perftools_LDADD) $(LIBUTILS)

succinct_trie_database_main_SOURCES = succinct_trie_database_main.cpp
succinct_trie_database_main_LDADD = $(boost_LDADD) $(LIBUTILS) $(perftools_LDADD) $(LIBUTILS)
