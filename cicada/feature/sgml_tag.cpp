//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/sgml_tag.hpp"
#include "cicada/parameter.hpp"
#include "cicada/weight_vector.hpp"

#include "utils/lexical_cast.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/compress_stream.hpp"
#include "utils/piece.hpp"
#include "utils/indexed_trie.hpp"

namespace cicada
{
  namespace feature
  {
    
    class SGMLTagImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;

      typedef cicada::HyperGraph hypergraph_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef attribute_set_type::attribute_type attribute_type;
      typedef feature_set_type::feature_type     feature_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef utils::indexed_trie<symbol_type, boost::hash<symbol_type>, std::equal_to<symbol_type>, std::allocator<symbol_type > > trie_type;
      typedef trie_type::id_type id_type;

      typedef std::vector<symbol_type, std::allocator<symbol_type> > stack_type;
      
      SGMLTagImpl() : feature("sgml-tag"), trie() {}
      
      inline
      bool is_sgml_tag(const symbol_type& symbol) const
      {
	return symbol.is_terminal() && symbol != vocab_type::EPSILON && symbol != vocab_type::BOS && symbol != vocab_type::EOS && symbol.is_sgml_tag();
      }
      
      void next_state(const symbol_type& tag, id_type& node, feature_set_type& features, int& penalty) const
      {
	trie_type& __trie = const_cast<trie_type&>(trie);
	
	if (tag.is_start_tag())
	  node = __trie.push(node, tag);
	else if (tag.is_end_tag()) {
	  if (node == __trie.root())
	    node = __trie.push(node, tag);
	  else {
	    const symbol_type& top = __trie[node];
	    
	    if (top.is_end_tag())
	      node = __trie.push(node, tag);
	    else {
	      penalty += 2 * (top.sgml_tag() != tag.sgml_tag());
	      node = __trie.pop(node);
	      
	      features[static_cast<const std::string&>(feature) + ':' + static_cast<const std::string&>(top) + '|' + static_cast<const std::string&>(tag)] += 1.0;
	    }
	  }
	}
      }
      
      void sgml_tag_score(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  const bool final) const
      {
	id_type node = trie.root();
	int penalty = 0;
	
	if (states.empty()) {
	  rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter)
	    if (is_sgml_tag(*riter))
	      next_state(*riter, node, features, penalty);
	} else {
	  stack_type& stack = const_cast<stack_type&>(stack_impl);
	  
	  int non_terminal_pos = 0;
	  rule_type::symbol_set_type::const_iterator riter_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator riter = edge.rule->rhs.begin(); riter != riter_end; ++ riter) {
	    if (riter->is_non_terminal()) {
	      const int __non_terminal_index = riter->non_terminal_index();
	      const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	      ++ non_terminal_pos;
	      
	      id_type node_antecedent = *reinterpret_cast<const id_type*>(states[antecedent_index]);
	      
	      if (node == trie.root())
		node = node_antecedent;
	      else if (node_antecedent != trie.root()) {
		stack.clear();
		
		while (node_antecedent != trie.root()) {
		  stack.push_back(trie[node_antecedent]);
		  node_antecedent = trie.pop(node_antecedent);
		}
		
		stack_type::const_reverse_iterator siter_end = stack.rend();
		for (stack_type::const_reverse_iterator siter = stack.rbegin(); siter != siter_end; ++ siter)
		  next_state(*siter, node, features, penalty);
	      }
	      
	    } else if (is_sgml_tag(*riter))
	      next_state(*riter, node, features, penalty);
	  }
	}

	if (final) {
	  id_type node_final = node;
	  while (node_final != trie.root()) {
	    features[static_cast<const std::string&>(feature) + ':' + static_cast<const std::string&>(trie[node_final])] += 1.0;
	    
	    ++ penalty;
	    node_final = trie.pop(node_final);
	  }
	}

	if (penalty == 0)
	  features.erase(feature);
	else
	  features[feature] = - penalty;
	
	*reinterpret_cast<id_type*>(state) = node;
      }
      
      feature_type feature;
      trie_type trie;
      stack_type stack_impl;
    };
    
    SGMLTag::SGMLTag(const std::string& parameter)
      : pimpl(new impl_type())
    {
      typedef boost::filesystem::path path_type;
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "sgml-tag")
	throw std::runtime_error("is this really sgmltag feature function? " + parameter);

      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	std::cerr << "WARNING: unsupported parameter for sgml-tag: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> sgml_tag(new impl_type());
      
      
      pimpl = sgml_tag.release();
      
      base_type::__state_size = sizeof(impl_type::id_type);
      base_type::__feature_name = "sgml-tag";
    }
    
    SGMLTag::~SGMLTag() { std::auto_ptr<impl_type> tmp(pimpl); }

    SGMLTag::SGMLTag(const SGMLTag& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    SGMLTag& SGMLTag::operator=(const SGMLTag& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void SGMLTag::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      pimpl->sgml_tag_score(state, states, edge, features, final);
    }
    
    void SGMLTag::apply_coarse(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void SGMLTag::apply_predict(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void SGMLTag::apply_scan(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     const int dot,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void SGMLTag::apply_complete(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {}
    
  };
};
