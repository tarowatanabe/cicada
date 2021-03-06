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
    class ComposeTreeCKY : public Operation
    {
    public:
      ComposeTreeCKY(const std::string& parameter,
		     const tree_grammar_type& __tree_grammar,
		     const grammar_type& __grammar,
		     const std::string& __goal,
		     const int __debug);
      
      void operator()(data_type& data) const;
      
      const tree_grammar_type& tree_grammar;
      const grammar_type&      grammar;
      tree_grammar_type tree_grammar_local;
      grammar_type      grammar_local;
      std::string goal;
      
      bool yield_source;
      bool frontier;
      bool unique_goal;
      
      int debug;
    };

    class ComposeTree : public Operation
    {
    public:
      ComposeTree(const std::string& parameter,
		  const tree_grammar_type& __tree_grammar,
		  const grammar_type& __grammar,
		  const std::string& __goal,
		  const int __debug);
      
      void operator()(data_type& data) const;
      
      const tree_grammar_type& tree_grammar;
      const grammar_type&      grammar;
      tree_grammar_type tree_grammar_local;
      grammar_type      grammar_local;
      std::string goal;
      
      bool yield_source;
      bool frontier;
      
      int debug;
    };

    class ComposeEarley : public Operation
    {
    public:
      ComposeEarley(const std::string& parameter,
		    const grammar_type& __grammar,
		    const std::string& __goal,
		    const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
      grammar_type grammar_local;
      std::string goal;
  
      bool yield_source;
      bool frontier;
  
      int debug;
    };

    class ComposeCKY : public Operation
    {
    public:
      ComposeCKY(const std::string& parameter,
		 const grammar_type& __grammar,
		 const std::string& __goal,
		 const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
      grammar_type grammar_local;
      std::string goal;
  
      bool yield_source;
      bool treebank;
      bool pos_mode;
      bool ordered;
      bool frontier;
      bool unique_goal;
  
      int debug;
    };

    class ComposeGrammar : public Operation
    {
    public:
      ComposeGrammar(const std::string& parameter,
		     const grammar_type& __grammar,
		     const std::string& __goal,
		     const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
      grammar_type grammar_local;
      std::string goal;

      bool yield_source;
      bool frontier;
      
      int debug;
    };


    class ComposePhrase : public Operation
    {
    public:
      ComposePhrase(const std::string& parameter,
		    const grammar_type& __grammar,
		    const std::string& __goal,
		    const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
      grammar_type grammar_local;
      std::string goal;
  
      int distortion;
      
      bool yield_source;
      bool frontier;
  
      int debug;
    };


    class ComposeAlignment : public Operation
    {
    public:
      ComposeAlignment(const std::string& parameter,
		       const grammar_type& __grammar,
		       const std::string& __goal,
		       const int __debug);
  
      void operator()(data_type& data) const;
  
      const grammar_type& grammar;
      grammar_type grammar_local;
      std::string goal;

      bool lattice_mode;
      bool forest_mode;
  
      int debug;
    };
    
    class ComposeDependency : public Operation
    {
    public:
      ComposeDependency(const std::string& parameter,
			const int __debug);
      
      void operator()(data_type& data) const;
      
      bool arc_standard;
      bool arc_eager;
      bool hybrid;
      bool top_down;
      bool degree2;
      
      int debug;
    };
  };
};


#endif
