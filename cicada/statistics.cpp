
#include "statistics.hpp"

namespace cicada
{
  
  std::ostream& operator<<(std::ostream& os, const Statistics::stat_type& stat)
  {
    os << "count: " << stat.count
       << " node: " << stat.node
       << " edge: " << stat.edge
       << " user-time: " << stat.user_time
       << " cpu-time: "  << stat.cpu_time;
    
    return os;
  }
  
  std::ostream& operator<<(std::ostream& os, const Statistics& stats)
  {
    
    
    return os;
  }
  
};
