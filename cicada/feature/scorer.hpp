// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__SCORER__HPP__
#define __CICADA__FEATURE__SCORER__HPP__ 1

#include <string>

#include <cicada/feature_function.hpp>
#include <cicada/eval.hpp>

// purely virtual scorer function...

namespace cicada
{
  namespace feature
  {
    class Scorer : public FeatureFunction
    {
    public:
      typedef cicada::eval::Score score_type;
      typedef score_type::score_ptr_type score_ptr_type;
      
      // this is scorer specific assign
      virtual void assign(const size_type& id,
			  const hypergraph_type& hypergraph,
			  const lattice_type& lattice,
			  const span_set_type& spans,
			  const sentence_set_type& targets,
			  const ngram_count_set_type& ngram_counts) = 0;

      virtual void assign(const score_ptr_type& score) = 0;
    };
    
  };
};

#endif
