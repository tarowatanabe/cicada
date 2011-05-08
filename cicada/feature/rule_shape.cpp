//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/rule_shape.hpp"
#include "cicada/parameter.hpp"

#include "utils/compact_trie_dense.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"
#include "utils/bithack.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    class RuleShapeImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;

      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef feature_function_type::rule_type rule_type;
      
      typedef utils::compact_trie_dense<int, feature_type, utils::hashmurmur<size_t>, std::equal_to<int>,
					std::allocator<std::pair<const int, feature_type> > > trie_type;
    
      typedef trie_type::id_type id_type;
    
      RuleShapeImpl() : trie(-1) {}
    
      void rule_shape_score(const edge_type& edge,
			    feature_set_type& features) const
      {
	id_type node = trie.root();
	
	int non_terminal_pos = 0;
	rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) {
	  int index = 0;
	  if (riter->is_non_terminal()) {
	    const int __non_terminal_index = riter->non_terminal_index();
	    const int non_terminal_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	    ++ non_terminal_pos;
	    
	    index = non_terminal_index + 1;
	  } 
	  
	  node = const_cast<trie_type&>(trie).insert(node, index);
	}
	
	if (trie[node].empty()) {
	  std::string feature("rule-shape");
	  
	  int non_terminal_pos = 0;
	  rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) {
	    int index = 0;
	    if (riter->is_non_terminal()) {
	      const int __non_terminal_index = riter->non_terminal_index();
	      const int non_terminal_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	      ++ non_terminal_pos;
	      
	      index = non_terminal_index + 1;
	    } 
	    
	    feature += ':' + utils::lexical_cast<std::string>(index);
	  }

	  const_cast<trie_type&>(trie).operator[](node) = feature;
	}
	
	features[trie[node]] = 1.0;
      }

      trie_type trie;
    };
    
    
    RuleShape::RuleShape(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "rule-shape")
	throw std::runtime_error("is this really rule-shape feature function? " + parameter);
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for rule-shape: " << piter->first << "=" << piter->second << std::endl;
      
      std::auto_ptr<impl_type> rule_shape_impl(new impl_type());
      
      base_type::__state_size = 0;
      base_type::__feature_name = "rule-shape";
      
      pimpl = rule_shape_impl.release();
    }
    
    RuleShape::~RuleShape() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    RuleShape::RuleShape(const RuleShape& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    { }
    
    RuleShape& RuleShape::operator=(const RuleShape& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void RuleShape::apply(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const
    {
      pimpl->rule_shape_score(edge, features);
    }

    void RuleShape::apply_coarse(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void RuleShape::apply_predict(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void RuleShape::apply_scan(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       const int dot,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {}
    
    void RuleShape::apply_complete(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {}
  };
};
