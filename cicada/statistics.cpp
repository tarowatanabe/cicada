
#include "statistics.hpp"

namespace cicada
{
  
  std::ostream& operator<<(std::ostream& os, const Statistics::stat_type& stat)
  {
    os << "count: " << stat.count
       << " node: " << stat.node
       << " edge: " << stat.edge
       << " user-time: " << stat.user_time
       << " cpu-time: "  << stat.cpu_time
       << " thread-time: " << stat.thread_time;
    
    return os;
  }
  
  std::ostream& operator<<(std::ostream& os, const Statistics& stats)
  {
    Statistics::const_iterator siter_end = stats.end();
    for (Statistics::const_iterator siter = stats.begin(); siter != siter_end; ++ siter)
      os << siter->first << ' ' << siter->second << '\n';
    
    return os;
  }
  
};
