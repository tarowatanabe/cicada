
#include <utils/config.hpp>
#include <utils/malloc_stats.hpp>

#if defined(HAVE_UNISTD_H)
#include <unistd.h>
#endif

#if defined(HAVE_SYS_RESOURCE_H)
#include <sys/resource.h>
#endif

#if defined(HAVE_MALLOC_MALLOC_H)
#include <malloc/malloc.h>
#endif

#if defined(HAVE_GOOGLE_MALLOC_EXTENSION_H) && defined(HAVE_GOOGLE_MALLOC_EXTENSION)
#include <google/malloc_extension.h>
#endif

#if defined(HAVE_JEMALLOC_H) && defined(HAVE_JEMALLOC_STATS)
#include <jemalloc.h>
#endif


namespace utils
{
  size_t malloc_stats::used()
  {
#if defined(HAVE_JEMALLOC_STATS)
    jemalloc_stats_t stats;
    jemalloc_stats(&stats);
    return stats.allocated;
#elif defined(HAVE_GOOGLE_MALLOC_EXTENSION)
    size_t num_used = 0;
    MallocExtension::instance()->GetNumericProperty("generic.current_allocated_bytes", &num_used);
    return num_used;
#elif defined(HAVE_MALLOC_MALLOC_H) && defined(HAVE_MALLOC_ZONE_STATISTICS)
    struct malloc_statistics_t stat;
    malloc_zone_statistics(NULL, &stat);
    return stat.size_in_use;
#elif defined(HAVE_SBRK)
    static char *memory_begin = reinterpret_cast<char*>(::sbrk(0));
    char *memory_end = reinterpret_cast<char*>(::sbrk(0));
    if (memory_end != ((char*)-1) && memofy_begin != ((char*)-1))
      return memory_end - memory_begin;
    else
      return 0;
#else
#warning "no reliable malloc statics..."
    return 0;
#endif
  }

  size_t malloc_stats::allocated()
  {
#if defined(HAVE_JEMALLOC_STATS)
    jemalloc_stats_t stats;
    jemalloc_stats(&stats);
    return stats.mapped;
#elif defined(HAVE_GOOGLE_MALLOC_EXTENSION)
    size_t num_allocated = 0;
    MallocExtension::instance()->GetNumericProperty("generic.heap_size", &num_allocated);
    return num_allocated;
#elif defined(HAVE_MALLOC_MALLOC_H) && defined(HAVE_MALLOC_ZONE_STATISTICS)
    struct malloc_statistics_t stat;
    malloc_zone_statistics(NULL, &stat);
    return stat.size_allocated;
#else
#warning "no reliable malloc statics..."
    return 0;
#endif
  }
};
