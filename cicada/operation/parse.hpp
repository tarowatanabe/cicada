// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__PARSE__HPP__
#define __CICADA__OPERATION__PARSE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class ParseCKY : public Operation
    {
    public:
      ParseCKY(const std::string& parameter,
	       const grammar_type& __grammar,
	       const std::string& __goal,
	       const int __debug);
  
      void operator()(data_type& data) const;
      
      void assign(const weight_set_type& __weights);
      
      const grammar_type& grammar;
  
      std::string goal;
  
      const weights_path_type* weights;
      const weight_set_type*   weights_assigned;
      int size;
      bool weights_one;
      
      bool yield_source;
      bool treebank;
      bool pos_mode;
      bool unique_goal;
      
      int debug;
    };
    
    class ParseAgenda : public Operation
    {
    public:
      ParseAgenda(const std::string& parameter,
		  const grammar_type& __grammar,
		  const std::string& __goal,
		  const int __debug);
      
      void operator()(data_type& data) const;
      
      void assign(const weight_set_type& __weights);
      
      const grammar_type& grammar;
      
      std::string goal;
      
      const weights_path_type* weights;
      const weight_set_type*   weights_assigned;
      int size;
      bool weights_one;
      
      bool yield_source;
      bool treebank;
      bool pos_mode;
      
      int debug;
    };

    class ParseCoarse : public Operation
    {
    private:
      typedef std::vector<grammar_type, std::allocator<grammar_type> > grammar_set_type;
      typedef std::vector<double, std::allocator<double> > threshold_set_type;
      
    public:
      ParseCoarse(const std::string& parameter,
		  const grammar_type& __grammar,
		  const std::string& __goal,
		  const int __debug);
      
      void operator()(data_type& data) const;
      
      void assign(const weight_set_type& __weights);
      
      grammar_set_type   grammars;
      threshold_set_type thresholds;
      
      const grammar_type& grammar;
      grammar_type grammar_local;
      std::string goal;
      
      const weights_path_type* weights;
      const weight_set_type*   weights_assigned;
      bool weights_one;
      
      bool yield_source;
      bool treebank;
      bool pos_mode;
      
      int debug;
    };

  };
};


#endif
