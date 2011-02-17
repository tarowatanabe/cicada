include mkinc_utf8.mak

noinst_LTLIBRARIES = libstemmer.la
noinst_HEADERS = $(snowball_headers)
libstemmer_la_SOURCES = $(snowball_sources)

noinst_PROGRAMS = stemwords

stemwords_SOURCES = examples/stemwords.c
stemwords_CPPFLAGS = -Iinclude -I. $(AM_CPPFLAGS) 
stemwords_LDADD = libstemmer.la
