// -*- mode: c++ -*-

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
    
    typedef operation_type::weight_set_type weight_set_type;
    typedef operation_type::grammar_type    grammar_type;
    typedef operation_type::model_type      model_type;
    typedef operation_type::path_type       path_type;
    
    typedef operation_type::hypergraph_type      hypergraph_type;
    typedef operation_type::lattice_type         lattice_type;
    typedef operation_type::span_set_type        span_set_type;
    typedef operation_type::sentence_type        sentence_type;
    typedef operation_type::sentence_set_type    sentence_set_type;
    typedef operation_type::ngram_count_set_type ngram_count_set_type;
    
    typedef boost::shared_ptr<operation_type> operation_ptr_type;
    typedef std::vector<operation_ptr_type, std::allocator<operation_ptr_type> > operation_ptr_set_type;
    
    typedef std::vector<std::string, std::allocator<std::string> > parameter_set_type;
    
  public:
    static std::string lists();
    
  public:
    template <typename Iterator>
    OperationSet(Iterator first, Iterator last,
		 const model_type& model,
		 const grammar_type& grammar,
		 const std::string& goal,
		 const std::string& non_terminal,
		 const bool insertion,
		 const bool deletion,
		 const bool input_id,
		 const bool input_lattice,
		 const bool input_forest,
		 const bool input_span,
		 const bool input_bitext,
		 const bool input_mpi,
		 const int debug)
    {
      parameter_set_type parameters(first, last);
      
      OperationSet(parameters, model, grammar, goal, non_terminal, insertion, deletion,
		   input_id, input_lattice, input_forest, input_span, input_bitext, input_mpi,
		   debug);
    }
    
    OperationSet(const parameter_set_type& parameters,
		 const model_type& model,
		 const grammar_type& grammar,
		 const std::string& goal,
		 const std::string& non_terminal,
		 const bool insertion,
		 const bool deletion,
		 const bool input_id,
		 const bool input_lattice,
		 const bool input_forest,
		 const bool input_span,
		 const bool input_bitext,
		 const bool input_mpi,
		 const int debug);
    
  public:
    void assign(const weight_set_type& weights);
    void operator()(const std::string& line);
    
    const output_data_type& get_output_data() const { return output_data; }
    const data_type& get_data() const { return data; }
    
  private:
    bool input_id;
    bool input_lattice;
    bool input_forest;
    bool input_span;
    bool input_bitext;
    bool input_mpi;
    
    output_data_type output_data;
    data_type        data;
    
    operation_ptr_set_type operations;
  };
};

#endif
