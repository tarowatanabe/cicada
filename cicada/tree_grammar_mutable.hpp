// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__TREE_GRAMMAR_MUTABLE__HPP__
#define __CICADA__TREE_GRAMMAR_MUTABLE__HPP__ 1

// mutable storage grammar...

#include <string>

#include <cicada/tree_transducer.hpp>

#include <boost/filesystem/path.hpp>

namespace cicada
{
  
  class TreeGrammarMutableImpl;

  class TreeGrammarMutable : public TreeTransducer
  {
  private:
    typedef TreeGrammarMutableImpl impl_type;
    
    typedef boost::filesystem::path path_type;
    
  public:
    // mutable grammar, use to encode rule-table generated by nicttm (probably...)
    // 
    //  [x]([x](a) b [x,1] c) ||| [x] (A [x,1] B) ||| score1 score2
    //
    // The parameter must be of form:
    // parameter = file-name:[key=value delimited by ',']*
    // where:
    // file-name: file name for the grammar, "-" for stdin
    // key,value: key, value pair... valid pairs are:
    //
    //    max-span = 15 : maximum non-terminals span
    // 
    //    feature0 = feature-name0
    //    feature1 = feature-name1
    //
    // if not supplied we will use rule-table-0, rule-table-1 etc.
    
    TreeGrammarMutable(const std::string& parameter);
    ~TreeGrammarMutable();
    
    TreeGrammarMutable(const TreeGrammarMutable& x);
    TreeGrammarMutable& operator=(const TreeGrammarMutable& x);

  private:
    TreeGrammarMutable() {}
    
  public:
    // virtual members
    transducer_ptr_type clone() const;
    
    edge_id_type edge(const symbol_type& symbol) const;
    edge_id_type edge(const symbol_set_type& symbols) const;
    edge_id_type edge(const symbol_type* first, const symbol_type* last) const;
    
    id_type root() const;
    id_type next(const id_type& node, const edge_id_type& edge) const;
    bool has_next(const id_type& node) const;
    const rule_pair_set_type& rules(const id_type& node) const;

    // grammar_mutable specific members
    void read(const std::string& parameter);
    void clear();
    

    void insert(const std::string& pattern);
    void insert(const rule_pair_type& rule_pair);
    
    
    void insert(const rule_type& source, const rule_type& target)
    {
      insert(rule_pair_type(rule_type::create(rule_type(source)), rule_type::create(rule_type(target))));
    }
    void insert(const rule_type& source, const rule_type& target, const feature_set_type& features)
    {
      insert(rule_pair_type(rule_type::create(rule_type(source)), rule_type::create(rule_type(target)), features));
    }
    void insert(const rule_type& source, const rule_type& target, const feature_set_type& features, const attribute_set_type& attributes)
    {
      insert(rule_pair_type(rule_type::create(rule_type(source)), rule_type::create(rule_type(target)), features, attributes));
    }
    void insert(const rule_ptr_type& source, const rule_ptr_type& target)
    {
      insert(rule_pair_type(source, target));
    }
    void insert(const rule_ptr_type& source, const rule_ptr_type& target, const feature_set_type& features)
    {
      insert(rule_pair_type(source, target, features));
    }
    void insert(const rule_ptr_type& source, const rule_ptr_type& target, const feature_set_type& features, const attribute_set_type& attributes)
    {
      insert(rule_pair_type(source, target, features, attributes));
    }
    
  private:
    impl_type* pimpl;
  };
  
};

#endif
