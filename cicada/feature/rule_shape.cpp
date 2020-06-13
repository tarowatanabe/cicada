//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/rule_shape.hpp"
#include "cicada/feature/feature_builder.hpp"

#include "cicada/parameter.hpp"

#include "utils/trie_compact.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"
#include "utils/bithack.hpp"
#include "utils/compact_map.hpp"
#include "utils/hashmurmur3.hpp"

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
      
      template <typename Tp>
      struct unassigned_key
      {
	Tp operator()() const { return Tp(-1); }
      };
      
      typedef utils::trie_compact<int, std::string,
				  unassigned_key<int>, 
				  boost::hash<int>, std::equal_to<int>,
				  std::allocator<std::pair<const int, std::string> > > trie_type;
      
      typedef std::vector<feature_type, std::allocator<feature_type> > cache_unigram_type;
      typedef std::vector<bool, std::allocator<bool> > checked_unigram_type;
      
      typedef trie_type::id_type id_type;
      
      typedef utils::compact_map<id_type, feature_type,
				 unassigned_key<id_type>, unassigned_key<id_type>,
				 boost::hash<id_type>, std::equal_to<id_type>,
				 std::allocator<std::pair<const id_type, feature_type> > > node_map_type;
      typedef std::deque<node_map_type, std::allocator<node_map_type> > cache_bigram_type;

      typedef FeatureBuilder feature_builder_type;
    
      RuleShapeImpl() : trie(), forced_feature(false) {}
    
      void rule_shape_score(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features)
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

	// root!
	if (node == trie.root())
	  node = const_cast<trie_type&>(trie).insert(node, -1);
	
	if (trie[node].empty()) {
	  std::string feature;
	  
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
	    
	    if (feature.empty())
	      feature = utils::lexical_cast<std::string>(index);
	    else
	      feature += '|' + utils::lexical_cast<std::string>(index);
	  }
	  
	  if (feature.empty())
	    feature = "<epsilon>";
	  
	  const_cast<trie_type&>(trie).operator[](node) = feature;
	}
	
	if (node >= cache_unigram.size())
	  cache_unigram.resize(node + 1);
	if (node >= cache_bigram.size())
	  cache_bigram.resize(node + 1);
	if (node >= checked_unigram.size())
	  checked_unigram.resize(node + 1, false);
	
	if (! checked_unigram[node]) {
	  checked_unigram[node] = true;

	  feature_builder.clear();
	  feature_builder << "rule-shape:" << trie[node];
	  
	  if (forced_feature || feature_builder.exists())
	    cache_unigram[node] = feature_builder;
	}
	
	if (! cache_unigram[node].empty())
	  features[cache_unigram[node]] += 1.0;
	
	*reinterpret_cast<id_type*>(state) = node;
	
	if (! states.empty()) {
	  const id_type& node_parent = node;
	  
	  int non_terminal_pos = 0;
	  rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
	  for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer) {
	    if (titer->is_non_terminal()) {
	      const int __non_terminal_index = titer->non_terminal_index();
	      const int antecedent_index = utils::bithack::branch(__non_terminal_index <= 0, non_terminal_pos, __non_terminal_index - 1);
	      ++ non_terminal_pos;
	      
	      const id_type& node_child = *reinterpret_cast<const id_type*>(states[antecedent_index]);
	      
	      std::pair<node_map_type::iterator, bool> result = cache_bigram[node_parent].insert(std::make_pair(node_child, feature_type()));
	      if (result.second) {
		feature_builder.clear();
		feature_builder << "rule-shape2:" << trie[node_parent] << "+" << trie[node_child];
		
		if (forced_feature || feature_builder.exists())
		  result.first->second = feature_builder;
	      }
	      
	      if (! result.first->second.empty())
		features[result.first->second] += 1.0;
	    }
	  }
	}
      }
      
      void clear()
      {
	checked_unigram.clear();
	cache_unigram.clear();
	cache_bigram.clear();
      }
      
      trie_type trie;
      
      checked_unigram_type checked_unigram;
      cache_unigram_type   cache_unigram;
      cache_bigram_type    cache_bigram;

      feature_builder_type feature_builder;
      
      bool forced_feature;
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
      
      std::unique_ptr<impl_type> rule_shape_impl(new impl_type());
      
      base_type::__state_size = sizeof(impl_type::id_type);
      base_type::__feature_name = "rule-shape";
      base_type::__sparse_feature = true;
      
      pimpl = rule_shape_impl.release();
    }
    
    RuleShape::~RuleShape() { std::unique_ptr<impl_type> tmp(pimpl); }
    
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
			  const bool final) const
    {
      feature_set_type feats;
 
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->rule_shape_score(state, states, edge, feats);
      
      features.update(feats, static_cast<const std::string&>(base_type::feature_name()));
    }

    void RuleShape::apply_coarse(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    void RuleShape::apply_predict(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  const bool final) const
    {
      apply(state, states, edge, features, final);
    }
    
    void RuleShape::apply_scan(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       const int dot,
			       feature_set_type& features,
			       const bool final) const
    {}
    
    void RuleShape::apply_complete(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   const bool final) const
    {}
    
    void RuleShape::initialize()
    {
      pimpl->clear();
    }
    
  };
};
