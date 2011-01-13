// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__OPERATION__COMPOSE__HPP__
#define __CICADA__OPERATION__COMPOSE__HPP__ 1

#include <iostream>

#include <cicada/operation.hpp>

namespace cicada
{
  namespace operation
  {
    class ComposeTree : public Operation
    {
    public:
      ComposeTree(const std::string& parameter,
		  const tree_grammar_type& __tree_grammar,
		  const grammar_type& __grammar,
		  const std::string& __goal,
		  const std::string& __non_terminal,
		  const bool __insertion,
		  const bool __deletion,
		  const bool __fallback,
		  const int __debug);
      
      void operator()(data_type& data) const;
      
      const tree_grammar_type& tree_grammar;
      const grammar_type&      grammar;
      
      std::string goal;
      std::string non_terminal;
      
      bool insertion;
      bool deletion;
      bool fallback;
      
      bool yield_source;
      
      int debug;
    };

    class ComposeEarley : public Operation
    {
    public:
      ComposeEarley(const std::string& parameter,
		    const grammar_type& __grammar,
		    const std::string& __goal,
		    const std::string& __non_terminal,
		    const bool __insertion,
		    const bool __deletion,
		    const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;
  
      bool insertion;
      bool deletion;

      bool yield_source;
  
      int debug;
    };

    class ComposeCKY : public Operation
    {
    public:
      ComposeCKY(const std::string& parameter,
		 const grammar_type& __grammar,
		 const std::string& __goal,
		 const std::string& __non_terminal,
		 const bool __insertion,
		 const bool __deletion,
		 const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;
  
      bool insertion;
      bool deletion;

      bool yield_source;
  
      int debug;
    };


    class ComposePhrase : public Operation
    {
    public:
      ComposePhrase(const std::string& parameter,
		    const grammar_type& __grammar,
		    const std::string& __goal,
		    const std::string& __non_terminal,
		    const bool __insertion,
		    const bool __deletion,
		    const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;
  
      bool insertion;
      bool deletion;

      int distortion;
      
      bool yield_source;
  
      int debug;
    };


    class ComposeAlignment : public Operation
    {
    public:
      ComposeAlignment(const std::string& parameter,
		       const grammar_type& __grammar,
		       const std::string& __goal,
		       const std::string& __non_terminal,
		       const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
  
      std::string goal;
      std::string non_terminal;

      bool lattice_mode;
      bool forest_mode;
  
      int debug;
    };
    
  };
};


#endif
