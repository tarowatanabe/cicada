// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION_SET__HPP__
#define __CICADA__OPERATION_SET__HPP__ 1

#include <string>
#include <vector>

#include <cicada/operation.hpp>


namespace cicada
{
  class OperationSet
  {
  public:
    typedef Operation     operation_type;
    
    typedef operation_type::data_type        data_type;
    typedef operation_type::output_data_type output_data_type;
    
    typedef operation_type::weight_set_type   weight_set_type;
    typedef operation_type::grammar_type      grammar_type;
    typedef operation_type::tree_grammar_type tree_grammar_type;
    typedef operation_type::model_type        model_type;
    
    typedef operation_type::path_type         path_type;
    
    typedef operation_type::hypergraph_type      hypergraph_type;
    typedef operation_type::lattice_type         lattice_type;
    typedef operation_type::span_set_type        span_set_type;
    typedef operation_type::alignment_type       alignment_type;
    typedef operation_type::sentence_type        sentence_type;
    typedef operation_type::sentence_set_type    sentence_set_type;
    typedef operation_type::ngram_count_set_type ngram_count_set_type;

    typedef operation_type::attribute_type  attribute_type;
    typedef operation_type::statistics_type statistics_type;
    
    typedef boost::shared_ptr<operation_type> operation_ptr_type;
    typedef std::vector<operation_ptr_type, std::allocator<operation_ptr_type> > operation_ptr_set_type;
    
    typedef std::vector<std::string, std::allocator<std::string> > parameter_set_type;
    
  public:
    static const char* lists();
    
  public:
    template <typename Iterator>
    OperationSet(Iterator first, Iterator last,
		 const model_type& model,
		 const grammar_type& grammar,
		 const tree_grammar_type& tree_grammar,
		 const std::string& goal,
		 const bool __input_id,
		 const bool __input_sentence,
		 const bool __input_lattice,
		 const bool __input_forest,
		 const bool __input_span,
		 const bool __input_alignment,
		 const bool __input_bitext,
		 const bool __input_mpi,
		 const int debug)
    {
      initialize(parameter_set_type(first, last), model, grammar, tree_grammar, goal, 
		 __input_id, __input_sentence, __input_lattice, __input_forest, __input_span, __input_alignment, __input_bitext, __input_mpi,
		 debug);
    }
    
    OperationSet(const parameter_set_type& parameters,
		 const model_type& model,
		 const grammar_type& grammar,
		 const tree_grammar_type& tree_grammar,
		 const std::string& goal,
		 const bool __input_id,
		 const bool __input_sentence,
		 const bool __input_lattice,
		 const bool __input_forest,
		 const bool __input_span,
		 const bool __input_alignment,
		 const bool __input_bitext,
		 const bool __input_mpi,
		 const int debug)
    {
      initialize(parameters, model, grammar, tree_grammar, goal, 
		 __input_id, __input_sentence, __input_lattice, __input_forest, __input_span, __input_alignment, __input_bitext, __input_mpi,
		 debug);
    }
    
  public:
    void clear();
    void assign(const weight_set_type& weights);
    
    void operator()(const std::string& line);
    
    const output_data_type& get_output_data() const { return output_data; }
    const data_type& get_data() const { return data; }

    const statistics_type& get_statistics() const { return statistics; }

    OperationSet clone() const
    {
      OperationSet __operations(*this);
      __operations.operations.clear();
      
      operation_ptr_set_type::const_iterator oiter_end = operations.end();
      for (operation_ptr_set_type::const_iterator oiter = operations.begin(); oiter != oiter_end; ++ oiter)
	__operations.operations.push_back(oiter->clone());
      
      return __operations;
    }

  private:
    void initialize(const parameter_set_type& parameters,
		    const model_type& model,
		    const grammar_type& grammar,
		    const tree_grammar_type& tree_grammar,
		    const std::string& goal,
		    const bool __input_id,
		    const bool __input_sentence,
		    const bool __input_lattice,
		    const bool __input_forest,
		    const bool __input_span,
		    const bool __input_alignment,
		    const bool __input_bitext,
		    const bool __input_mpi,
		    const int debug);
    
  private:
    bool input_id;
    bool input_sentence;
    bool input_lattice;
    bool input_forest;
    bool input_span;
    bool input_alignment;
    bool input_bitext;
    bool input_mpi;
    
    output_data_type output_data;
    data_type        data;
    
    operation_ptr_set_type operations;
    statistics_type        statistics;
  };
};

#endif
