
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
chinese_restaurant_process.hpp \
chunk_vector.hpp \
compact_trie.hpp \
compact_trie_dense.hpp \
compact_trie_set.hpp \
compact_trie_dense_set.hpp \
compress_stream.hpp \
config.hpp \
config.h \
dense_hashtable.hpp \
dense_hash_map.hpp \
dense_hash_set.hpp \
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
mulvector2.hpp \
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
sampler.hpp \
simple_vector.hpp \
slice_sampler.hpp \
small_vector.hpp \
space_separator.hpp \
sparse_hash_map.hpp \
sparse_hash_set.hpp \
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
unordered_map.hpp \
unordered_set.hpp \
vector_map.hpp \
vector_set.hpp \
vector2.hpp \
vector3.hpp \
vector2_aligned.hpp \
vector3_aligned.hpp \
vertical_coded_device.hpp \
vertical_coded_vector.hpp

noinst_PROGRAMS = \
arc_list_main \
b_heap_main \
base64_main \
bithack_main \
bitpack_main \
bit_vector_main \
byte_aligned_code_main \
byte_aligned_vector_main \
chinese_restaurant_process_main \
chunk_vector_main \
compress_stream_main \
dense_map_main \
double_base64_main \
filesystem_main \
group_aligned_code_main \
icu_filter_main \
indexed_set_main \
istream_line_iterator_main \
lexical_cast_main \
lockfree_list_queue_main \
lockfree_queue_main \
malloc_stats_main \
map_file_allocator_main \
packed_vector_main \
piece_main \
program_options_main \
random_seed_main \
simple_vector_main \
small_vector_main \
sparse_map_main \
static_allocator_main \
subprocess_main \
succinct_vector_main \
tempfile_main \
json_string_main \
vector_map_main \
vector_set_main \
vector2_main \
vertical_coded_vector_main

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

group_aligned_code_main.cpp : group_aligned_code.hpp group_aligned_code_offsets.hpp

group_aligned_code.hpp : group_aligned_code_offsets.hpp

group_aligned_code_offsets.hpp : group_aligned_code_offsets.py
	./group_aligned_code_offsets.py > group_aligned_code_offsets.hpp

arc_list_main_SOURCES = arc_list_main.cpp
arc_list_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS) 
arc_list_main_LDADD = $(BOOST_THREAD_LIBS)

b_heap_main_SOURCES = b_heap_main.cpp
b_heap_main_LDADD = $(LIBUTILS)

base64_main_SOURCES = base64_main.cpp
#base64_main_LDADD = $(LIBUTILS)

bithack_main_SOURCES = bithack_main.cpp

bitpack_main_SOURCES = bitpack_main.cpp
bitpack_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
bitpack_main_LDADD = $(LIBUTILS) $(BOOST_THREAD_LIBS)

bit_vector_main_SOURCES = bit_vector_main.cpp

byte_aligned_code_main_SOURCES = byte_aligned_code_main.cpp

byte_aligned_vector_main_SOURCES = byte_aligned_vector_main.cpp

chinese_restaurant_process_main_SOURCES = chinese_restaurant_process_main.cpp

chunk_vector_main_SOURCES = chunk_vector_main.cpp
chunk_vector_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
chunk_vector_main_LDADD = $(BOOST_THREAD_LIBS)

compress_stream_main_SOURCES = compress_stream_main.cpp
compress_stream_main_LDADD = $(LIBUTILS)

dense_map_main_SOURCES = dense_map_main.cpp
dense_map_main_LDADD = $(LIBUTILS)

double_base64_main_SOURCES = double_base64_main.cpp
double_base64_main_LDADD = $(LIBUTILS)

filesystem_main_SOURCES = filesystem_main.cpp
filesystem_main_LDADD = $(LIBUTILS)

group_aligned_code_main_SOURCES = group_aligned_code_main.cpp

icu_filter_main_SOURCES = icu_filter_main.cpp
icu_filter_main_LDADD = $(LIBUTILS) $(BOOST_IOSTREAMS_LDFLAGS)

indexed_set_main_SOURCES = indexed_set_main.cpp
indexed_set_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS) $(BOOST_THREAD_LIBS)

istream_line_iterator_main_SOURCES = istream_line_iterator_main.cpp

lexical_cast_main_SOURCES = lexical_cast_main.cpp
lexical_cast_main_LDADD = $(LIBUTILS)

lockfree_list_queue_main_SOURCES = lockfree_list_queue_main.cpp
lockfree_list_queue_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
lockfree_list_queue_main_LDADD = $(BOOST_THREAD_LIBS)

lockfree_queue_main_SOURCES = lockfree_queue_main.cpp
lockfree_queue_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
lockfree_queue_main_LDADD = $(BOOST_THREAD_LIBS)

malloc_stats_main_SOURCES = malloc_stats_main.cpp
malloc_stats_main_LDADD = $(LIBUTILS)

map_file_allocator_main_SOURCES = map_file_allocator_main.cpp
map_file_allocator_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
map_file_allocator_main_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_THREAD_LIBS) $(LIBUTILS) $(perftools_LDADD)

packed_vector_main_SOURCES = packed_vector_main.cpp
packed_vector_main_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
packed_vector_main_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_THREAD_LIBS) $(LIBUTILS)

piece_main_SOURCES = piece_main.cpp
piece_main_LDADD = $(LIBUTILS)

program_options_main_SOURCES = program_options_main.cpp
program_options_main_LDADD = $(BOOST_PROGRAM_OPTIONS_LDFLAGS) $(BOOST_PROGRAM_OPTIONS_LIBS)

simple_vector_main_SOURCES = simple_vector_main.cpp
simple_vector_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
simple_vector_main_LDADD = $(BOOST_THREAD_LIBS)

small_vector_main_SOURCES = small_vector_main.cpp
small_vector_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
small_vector_main_LDADD = $(BOOST_THREAD_LIBS)

static_allocator_main_SOURCES = static_allocator_main.cpp
static_allocator_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
static_allocator_main_LDADD = $(BOOST_THREAD_LIBS)

sparse_map_main_SOURCES = sparse_map_main.cpp
sparse_map_main_LDADD = $(LIBUTILS)

random_seed_main_SOURCES = random_seed_main.cpp
random_seed_main_LDADD = $(LIBUTILS)

subprocess_main_SOURCES = subprocess_main.cpp
subprocess_main_LDFLAGS = $(BOOST_IOSTREAMS_LDFLAGS) $(BOOST_THREAD_LDFLAGS) $(BOOST_FILESYSTEM_LDFLAGS)
subprocess_main_LDADD = $(BOOST_IOSTREAMS_LIBS) $(BOOST_FILESYSTEM_LIBS) $(BOOST_THREAD_LIBS)

succinct_vector_main_SOURCES = succinct_vector_main.cpp
succinct_vector_main_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
succinct_vector_main_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_THREAD_LIBS)

tempfile_main_SOURCES = tempfile_main.cpp
tempfile_main_LDFLAGS = $(BOOST_THREAD_LDFLAGS)
tempfile_main_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_THREAD_LIBS) $(LIBUTILS)

json_string_main_SOURCES = json_string_main.cpp

vector2_main_SOURCES = vector2_main.cpp

vector_map_main_SOURCES = vector_map_main.cpp

vector_set_main_SOURCES = vector_set_main.cpp

vertical_coded_vector_main_SOURCES = vertical_coded_vector_main.cpp
vertical_coded_vector_main_LDFLAGS = $(BOOST_FILESYSTEM_LDFLAGS) $(BOOST_THREAD_LDFLAGS)
vertical_coded_vector_main_LDADD = $(BOOST_FILESYSTEM_LIBS) $(BOOST_IOSTREAMS_LIBS) $(BOOST_THREAD_LIBS)

