// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__HPP__
#define __CICADA__OPERATION__HPP__ 1

#include <iostream>
#include <string>

#include <cicada/hypergraph.hpp>
#include <cicada/lattice.hpp>
#include <cicada/sentence.hpp>
#include <cicada/sentence_vector.hpp>
#include <cicada/alignment.hpp>
#include <cicada/dependency.hpp>
#include <cicada/span_vector.hpp>
#include <cicada/ngram_count_set.hpp>
#include <cicada/span_vector.hpp>
#include <cicada/weight_vector.hpp>
#include <cicada/grammar.hpp>
#include <cicada/tree_grammar.hpp>
#include <cicada/model.hpp>
#include <cicada/vocab.hpp>
#include <cicada/attribute.hpp>
#include <cicada/statistics.hpp>

#include <boost/filesystem/path.hpp>
#include <boost/shared_ptr.hpp>

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
    typedef cicada::Alignment      alignment_type;
    typedef cicada::Dependency     dependency_type;
    typedef cicada::Sentence       sentence_type;
    typedef cicada::SentenceVector sentence_set_type;
    typedef cicada::NGramCountSet  ngram_count_set_type;
    
    typedef cicada::Symbol symbol_type;
    typedef cicada::Symbol word_type;
    typedef cicada::Vocab  vocab_type;

    typedef cicada::Grammar         grammar_type;
    typedef cicada::TreeGrammar     tree_grammar_type;
    typedef cicada::Model           model_type;
    typedef cicada::FeatureFunction feature_function_type;

    typedef cicada::Attribute  attribute_type;
    typedef cicada::Statistics statistics_type;
    
    typedef feature_function_type::feature_function_ptr_type feature_function_ptr_type;
    
    typedef boost::filesystem::path path_type;    
    
    typedef Operation base_type;
    typedef Operation operation_type;
    typedef boost::shared_ptr<operation_type> operation_ptr_type;

    struct weights_path_type
    {
      weight_set_type weights;
      path_type       path;
      
      weights_path_type() : weights(), path() {}
      weights_path_type(const path_type& __path) : weights(), path(__path) {}
    };
    
    struct Data
    {
      size_type id;
      
      hypergraph_type      hypergraph;
      lattice_type         lattice;
      span_set_type        spans;
      alignment_type       alignment;
      dependency_type      dependency;
      sentence_set_type    targets;
      ngram_count_set_type ngram_counts;
      statistics_type      statistics;
    };

    struct OutputData
    {
      boost::shared_ptr<std::ostream> os;
      std::string                     buffer;
      
      path_type file;
      path_type directory;
      
      bool use_buffer;
    };
    
    typedef Data       data_type;
    typedef OutputData output_data_type;
    
    Operation() : name() {}
    Operation(const attribute_type& __name) : name(__name) {}
    virtual ~Operation() {}
    
    virtual operation_ptr_type clone() const { return operation_ptr_type(); }
    virtual void operator()(data_type& data) const = 0;
    virtual void assign(const weight_set_type& weights) {}
    virtual void clear() {};
    
    static const weights_path_type& weights();
    static const weights_path_type& weights(const path_type& path);
    
  protected:
    attribute_type name;
  };
};

#endif
