// -*- mode: c++ -*-

#ifndef __CICADA__OPERATION__HPP__
#define __CICADA__OPERATION__HPP__ 1

#include <cicada/hypergraph.hpp>
#include <cicada/lattice.hpp>
#include <cicada/sentence.hpp>
#include <cicada/sentence_vector.hpp>
#include <cicada/span_vector.hpp>
#include <cicada/ngram_count_set.hpp>
#include <cicada/span_vector.hpp>
#include <cicada/weight_vector.hpp>

#include <boost/filesystem/path.hpp>

namespace cicada
{

  struct Operation
  {
    typedef size_t    size_type;
    typedef ptrdiff_t difference_type;
    
    typedef cicada::WeightVector<double>   weight_set_type;
    
    typedef cicada::HyperGraph     hypergraph_type;
    typedef cicada::Lattice        lattice_type;
    typedef cicada::SpanVector     span_set_type;
    typedef cicada::Sentence       sentence_type;
    typedef cicada::SentenceVector sentence_set_type;
    typedef cicada::NGramCountSet  ngram_count_set_type;
    
    typedef boost::filesystem::path path_type;    
    
    typedef Operation base_type;
    
    struct Data
    {
      size_type id;
      
      hypergraph_type      hypergraph;
      lattice_type         lattice;
      span_set_type        spans;
      sentence_set_type    targets;
      ngram_count_set_type ngram_counts;
    };
    typedef Data data_type;
    
    
    Operation() {}
    virtual ~Operation() {}
    
    virtual void operator()(data_type& data) const = 0;
    virtual void assign(const weight_set_type& weights) {}
    virtual void clear() {};
    
    static const weight_set_type& weights(const path_type& path);
  };
  
  class OperationSet
  {
  public:
    typedef Operation     operation_type;
    
    typedef operation_type::data_type data_type;
    typedef operation_type::weight_set_type weight_set_type;
    
    typedef boost::shared_ptr<operation_type> operation_ptr_type;
    typedef std::vector<operation_ptr_type, std::allocator<operation_ptr_type> > operation_ptr_set_type;
    
  public:
    static std::string lists();

  public:
    template <typename Iterator>
    OperationSet(Iterator first, Iterator last) {}
    
  public:
    void assign(const weight_set_type& weights);
    void operator()(const std::string& line);
  };
  
};

#endif
