//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "line_search.hpp"

namespace cicada
{
  namespace optimize
  {
    const double LineSearch::interval_min = 1e-4;
    const double LineSearch::interval_offset_lower = 200.0;
    const double LineSearch::interval_offset_upper = 0.2;

    double LineSearch::value_min = - 100.0;
    double LineSearch::value_max =   100.0;
  }
};

