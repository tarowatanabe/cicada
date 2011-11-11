
AM_CPPFLAGS = -I$(top_builddir) $(BOOST_CPPFLAGS) $(ICU_CPPFLAGS) $(SPARSEHASH_CPPFLAGS) $(TCMALLOC_CPPFLAGS) $(JEMALLOC_CPPFLAGS) $(PTHREAD_FLAGS)

perftools_LDADD = $(JEMALLOC_LDFLAGS) $(TCMALLOC_LDFLAGS) $(PROFILER_LDFLAGS)

LIBUTILS=libutils.la

lib_LTLIBRARIES = libutils.la

libutils_la_SOURCES = \
	bitpack8.cpp \
	bitpack16.cpp \
	bitpack32.cpp \
	bitpack64.cpp \
	icu_filter.cpp \
	malloc_stats.cpp \
	map_file_allocator.cpp \
	tempfile.cpp

libutils_la_LDFLAGS = -version-info $(CICADA_LTVERSION)

libutils_la_LIBADD = \
	$(ICU_LDFLAGS) \
	$(BOOST_THREAD_LDFLAGS) $(BOOST_THREAD_LIBS) \
	$(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_IOSTREAMS_LIBS) \
	$(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_FILESYSTEM_LIBS) \
	$(BOOST_SYSTEM_LDFLAGS) $(BOOST_SYSTEM_LIBS) \
	$(JEMALLOC_LDFLAGS) \
	$(TCMALLOC_LDFLAGS)

libutilsincludedir = $(pkgincludedir)/utils
libutilsinclude_HEADERS = \
alloc_vector.hpp \
allocinfo_allocator.hpp \
arc_list.hpp \
array_power2.hpp \
async_device.hpp \
atomicop.hpp \
b_heap.hpp \
base64.hpp \
bithack.hpp \
bitpack.hpp \
bit_vector.hpp \
bounded_queue.hpp \
byte_aligned_code.hpp \
byte_aligned_vector.hpp \
byte_aligned_delta_vector.hpp \
byte_aligned_pair_delta_vector.hpp \
c_string_generator.hpp \
c_string_parser.hpp \
chart.hpp \
chunk_vector.hpp \
compact_trie.hpp \
compact_trie_dense.hpp \
compact_trie_set.hpp \
compact_trie_dense_set.hpp \
compress_stream.hpp \
config.hpp \
config.h \
dense_hashtable.hpp \
dense_map.hpp \
dense_set.hpp \
double_base64_generator.hpp \
double_base64_parser.hpp \
filesystem.hpp \
group_aligned_code.hpp \
group_aligned_code_offsets.hpp \
hashmurmur.hpp \
icu_filter.hpp \
indexed_hashtable.hpp \
indexed_map.hpp \
indexed_set.hpp \
indexed_trie.hpp \
istream_line_iterator.hpp \
lexical_cast.hpp \
lockfree_list_queue.hpp \
lockfree_queue.hpp \
malloc_stats.hpp \
map_file.hpp \
map_file_allocator.hpp \
mathop.hpp \
memory.hpp \
mpi.hpp \
mpi_allocator.hpp \
mpi_device.hpp \
mpi_device_bcast.hpp \
mpi_stream.hpp \
mpi_stream_simple.hpp \
mpi_traits.hpp \
packed_device.hpp \
packed_vector.hpp \
piece.hpp \
program_options.hpp \
python_string_parser.hpp \
python_string_generator.hpp \
random_seed.hpp \
regex_group.hpp \
repository.hpp \
resource.hpp \
sgi_hash_map.hpp \
sgi_hash_set.hpp \
simple_vector.hpp \
space_separator.hpp \
sparse_hashtable.hpp \
sparse_map.hpp \
sparse_set.hpp \
spinlock.hpp \
static_allocator.hpp \
std_heap.hpp \
subprocess.hpp \
succinct_vector.hpp \
symbol.hpp \
symbol_sequence.hpp \
symbol_tree.hpp \
tempfile.hpp \
thread_specific_ptr.hpp \
json_string_generator.hpp \
json_string_parser.hpp \
vector_map.hpp \
vector_set.hpp \
vector2.hpp \
vector3.hpp \
vector2_aligned.hpp \
vector3_aligned.hpp \
vertical_coded_device.hpp \
vertical_coded_vector.hpp

noinst_PROGRAMS = \
sample_arc_list \
sample_b_heap \
sample_base64 \
sample_bithack \
sample_bitpack \
sample_bit_vector \
sample_byte_aligned_code \
sample_byte_aligned_vector \
sample_chunk_vector \
sample_compress_stream \
sample_dense_map \
sample_double_base64 \
sample_filesystem \
sample_group_aligned_code \
sample_icu_filter \
sample_indexed_set \
sample_istream_line_iterator \
sample_lexical_cast \
sample_lockfree_list_queue \
sample_lockfree_queue \
sample_malloc_stats \
sample_map_file_allocator \
sample_packed_vector \
sample_piece \
sample_program_options \
sample_simple_vector \
sample_sparse_map \
sample_static_allocator \
sample_subprocess \
sample_succinct_vector \
sample_tempfile \
sample_json_string \
sample_vector_map \
sample_vector_set \
sample_vector2 \
sample_vertical_coded_vector

dist_noinst_SCRIPTS = \
bitpack_gen.py \
group_aligned_code_offsets.py

IGNORES = \
bitpack16_impl.hpp \
bitpack32_impl.hpp \
bitpack64_impl.hpp \
bitpack8_impl.hpp \
config.h \
config.h.in \
group_aligned_code_offsets.hpp \
stamp-h1

noinst_DATA = .gitignore

.gitignore: Makefile.am
	cat $(top_builddir)/.gitignore > .gitignore
	for bin in $(PROGRAMS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(SCRIPTS); do \
	  echo $$bin >> .gitignore; \
	done; \
	for bin in $(IGNORES); do \
	  echo $$bin >> .gitignore; \
	done

noinst_HEADERS = \
bitpack_impl.hpp

bitpack8.cpp  : bitpack_impl.hpp bitpack8_impl.hpp
bitpack16.cpp  : bitpack_impl.hpp bitpack16_impl.hpp
bitpack32.cpp : bitpack_impl.hpp bitpack32_impl.hpp
bitpack64.cpp : bitpack_impl.hpp bitpack64_impl.hpp

bitpack8_impl.hpp : bitpack_gen.py
	./bitpack_gen.py --byte-size 1 > bitpack8_impl.hpp

bitpack16_impl.hpp : bitpack_gen.py
	./bitpack_gen.py --byte-size 2 > bitpack16_impl.hpp

bitpack32_impl.hpp : bitpack_gen.py
	./bitpack_gen.py --byte-size 4 > bitpack32_impl.hpp

bitpack64_impl.hpp : bitpack_gen.py
	./bitpack_gen.py --byte-size 8 > bitpack64_impl.hpp

sample_group_aligned_code.cpp : group_aligned_code.hpp group_aligned_code_offsets.hpp

group_aligned_code.hpp : group_aligned_code_offsets.hpp

group_aligned_code_offsets.hpp : group_aligned_code_offsets.py
	./group_aligned_code_offsets.py > group_aligned_code_offsets.hpp

sample_arc_list_SOURCES = sample_arc_list.cpp
sample_arc_list_LDFLAGS = $(BOOST_THREAD_LDFLAGS) 
sample_arc_list_LDADD = $(BOOST_THREAD_LIBS)

sample_b_heap_SOURCES = sample_b_heap.cpp
sample_b_heap_LDADD = $(LIBUTILS)

sample_base64_SOURCES = sample_base64.cpp
#sample_base64_LDADD = $(LIBUTILS)

sample_bithack_SOURCES = sample_bithack.cpp

sample_bitpack_SOURCES = sample_bitpack.cpp
sample_bitpack_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
sample_bitpack_LDADD = $(LIBUTILS) $(BOOST_THREAD_LIBS)

sample_bit_vector_SOURCES = sample_bit_vector.cpp

sample_byte_aligned_code_SOURCES = sample_byte_aligned_code.cpp

sample_byte_aligned_vector_SOURCES = sample_byte_aligned_vector.cpp

sample_chunk_vector_SOURCES = sample_chunk_vector.cpp
sample_chunk_vector_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
sample_chunk_vector_LDADD = $(BOOST_THREAD_LIBS)

sample_compress_stream_SOURCES = sample_compress_stream.cpp
sample_compress_stream_LDADD = $(LIBUTILS)

sample_dense_map_SOURCES = sample_dense_map.cpp
sample_dense_map_LDADD = $(LIBUTILS)

sample_double_base64_SOURCES = sample_double_base64.cpp
sample_double_base64_LDADD = $(LIBUTILS)

sample_filesystem_SOURCES = sample_filesystem.cpp
sample_filesystem_LDADD = $(LIBUTILS)

sample_group_aligned_code_SOURCES = sample_group_aligned_code.cpp

sample_icu_filter_SOURCES = sample_icu_filter.cpp
sample_icu_filter_LDADD = $(LIBUTILS) $(BOOST_IOSTREAMS_LDFLAGS)

sample_indexed_set_SOURCES = sample_indexed_set.cpp
sample_indexed_set_LDFLAGS = $(BOOST_THREAD_LDFLAGS) $(BOOST_THREAD_LIBS)

sample_istream_line_iterator_SOURCES = sample_istream_line_iterator.cpp

sample_lexical_cast_SOURCES = sample_lexical_cast.cpp
sample_lexical_cast_LDADD = $(LIBUTILS)

sample_lockfree_list_queue_SOURCES = sample_lockfree_list_queue.cpp
sample_lockfree_list_queue_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
sample_lockfree_list_queue_LDADD = $(BOOST_THREAD_LIBS) $(BOOST_SYSTEM_LIBS)

sample_lockfree_queue_SOURCES = sample_lockfree_queue.cpp
sample_lockfree_queue_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
sample_lockfree_queue_LDADD = $(BOOST_THREAD_LIBS) $(BOOST_SYSTEM_LIBS)

sample_malloc_stats_SOURCES = sample_malloc_stats.cpp
sample_malloc_stats_LDADD = $(LIBUTILS)

sample_map_file_allocator_SOURCES = sample_map_file_allocator.cpp
sample_map_file_allocator_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
sample_map_file_allocator_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_THREAD_LIBS) $(LIBUTILS) $(perftools_LDADD)


sample_packed_vector_SOURCES = sample_packed_vector.cpp
sample_packed_vector_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
sample_packed_vector_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_THREAD_LIBS) $(LIBUTILS)

sample_piece_SOURCES = sample_piece.cpp
sample_piece_LDADD = $(LIBUTILS)

sample_program_options_SOURCES = sample_program_options.cpp
sample_program_options_LDADD = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

sample_simple_vector_SOURCES = sample_simple_vector.cpp
sample_simple_vector_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
sample_simple_vector_LDADD = $(BOOST_THREAD_LIBS)

sample_static_allocator_SOURCES = sample_static_allocator.cpp
sample_static_allocator_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
sample_static_allocator_LDADD = $(BOOST_THREAD_LIBS)

sample_sparse_map_SOURCES = sample_sparse_map.cpp
sample_sparse_map_LDADD = $(LIBUTILS)

sample_subprocess_SOURCES = sample_subprocess.cpp
sample_subprocess_LDFLAGS = $(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_THREAD_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS)
sample_subprocess_LDADD = $(BOOST_IOSTREAMS_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_THREAD_LIBS)

sample_succinct_vector_SOURCES = sample_succinct_vector.cpp
sample_succinct_vector_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
sample_succinct_vector_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_THREAD_LIBS)

sample_tempfile_SOURCES = sample_tempfile.cpp
sample_tempfile_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
sample_tempfile_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_THREAD_LIBS) $(BOOST_SYSTEM_LIBS) $(LIBUTILS)

sample_json_string_SOURCES = sample_json_string.cpp

sample_vector2_SOURCES = sample_vector2.cpp

sample_vector_map_SOURCES = sample_vector_map.cpp

sample_vector_set_SOURCES = sample_vector_set.cpp

sample_vertical_coded_vector_SOURCES = sample_vertical_coded_vector.cpp
sample_vertical_coded_vector_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
sample_vertical_coded_vector_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_SYSTEM_LIBS) $(BOOST_THREAD_LIBS)

