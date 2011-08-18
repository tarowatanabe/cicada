AM_CPPFLAGS = -I.. $(BOOST_CPPFLAGS)

LIBUTILS=$(top_builddir)/utils/libutils.la

perftools_LDADD = $(JEMALLOC_LDFLAGS) $(TCMALLOC_LDFLAGS) $(PROFILER_LDFLAGS)

boost_LDADD = \
	$(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS) \
	$(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_IOSTREAMS_LIBS) \
	$(BOOST_THREAD_LDFLAGS) $(BOOST_THREAD_LIBS) \
	$(BOOST_SYSTEM_LDFLAGS) $(BOOST_SYSTEM_LIBS)

succinct_dbincludedir = $(pkgincludedir)/succinct_db
succinct_dbinclude_HEADERS = \
	succinct_trie.hpp \
	succinct_hash.hpp \
	succinct_trie_db.hpp \
	succinct_trie_database.hpp \
	succinct_hash_db.hpp

noinst_PROGRAMS = \
	sample_succinct_trie \
	sample_succinct_hash \
	sample_succinct_trie_db \
	sample_succinct_trie_database

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done

sample_succinct_trie_SOURCES = sample_succinct_trie.cpp
sample_succinct_trie_LDADD = $(boost_LDADD) $(LIBUTILS) $(perftools_LDADD)

sample_succinct_hash_SOURCES = sample_succinct_hash.cpp
sample_succinct_hash_LDADD = $(boost_LDADD) $(LIBUTILS) $(perftools_LDADD)

sample_succinct_trie_db_SOURCES = sample_succinct_trie_db.cpp
sample_succinct_trie_db_LDADD = $(boost_LDADD) $(LIBUTILS) $(perftools_LDADD) $(LIBUTILS)

sample_succinct_trie_database_SOURCES = sample_succinct_trie_database.cpp
sample_succinct_trie_database_LDADD = $(boost_LDADD) $(LIBUTILS) $(perftools_LDADD) $(LIBUTILS)
