// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__HPP__
#define __CICADA__OPERATION__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/lattice.hpp>
#include <cicada/sentence.hpp>
#include <cicada/ngram_count_set.hpp>
#include <cicada/span_vector.hpp>

namespace cicada
{

  struct Operation
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::WeightVector<double>   weight_set_type;
    
    typedef cicada::HyperGraph    hypergraph_type;
    typedef cicada::Lattice       lattice_type;
    typedef cicada::Sentence      sentence_type;
    typedef cicada::SpanVector    span_set_type;
    typedef cicada::NGramCountSet ngram_count_set_type;
    typedef std::vector<sentence_type, std::allocator<sentence_type> > sentence_set_type;
    
    typedef Operation base_type;
    
    struct Parameter
    {
      size_type id;
      
      lattice_type         lattice;
      span_set_type        spans;
      sentence_set_type    targets;
      hypergraph_type      hypergraph;
      ngram_count_set_type ngram_counts;
    };
    typedef Parameter parameter_type;
    
    
    Operation() {}
    virtual ~Operation() {}
    
    virtual void operator()(data_type& data) const = 0;
    
    virtual void assign(const weight_set_type& weights) {}
    
    virtual void clear() {};
  };

  class OperationSet
  {
    
    
  };
  
};

#endif
