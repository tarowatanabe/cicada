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
	       const std::string& __non_terminal,
	       const bool __insertion,
	       const bool __deletion,
	       const int __debug);
  
      void operator()(data_type& data) const;
      
      void assign(const weight_set_type& __weights);
      
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;
  
      bool insertion;
      bool deletion;
      
      const weight_set_type* weights;
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
		  const std::string& __non_terminal,
		  const bool __insertion,
		  const bool __deletion,
		  const int __debug);
      
      void operator()(data_type& data) const;
      
      void assign(const weight_set_type& __weights);
      
      const grammar_type& grammar;
      
      std::string goal;
      std::string non_terminal;
      
      bool insertion;
      bool deletion;

      const weight_set_type* weights;
      int size;
      bool weights_one;
      
      bool yield_source;
      bool treebank;
      bool pos_mode;
      
      int debug;
    };

  };
};


#endif
