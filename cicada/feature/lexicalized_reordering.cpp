//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/lexicalized_reordering.hpp"
#include "cicada/parameter.hpp"

#include "utils/indexed_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/simple_vector.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    
    class LexicalizedReorderingImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      typedef cicada::Lattice  lattice_type;

      typedef size_t    size_type;
      typedef ptrdiff_t difference_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      
      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef std::vector<feature_type, std::allocator<feature_type> >     feature_list_type;
      typedef std::vector<attribute_type, std::allocator<attribute_type> > attribute_list_type;
      
      typedef utils::simple_vector<float, std::allocator<float> > feature_cache_type;
      
      struct feature_cache_hash_type : public utils::hashmurmur<size_t>
      {
	size_t operator()(const feature_cache_type& x) const
	{
	  return utils::hashmurmur<size_t>::operator()(x.begin(), x.end(), 0);
	}
      };
      typedef utils::indexed_set<feature_cache_type, feature_cache_hash_type, std::equal_to<feature_cache_type>,
				 std::allocator<feature_cache_type> > feature_cache_states_type;
      typedef feature_cache_states_type::index_type state_type;
      
      LexicalizedReorderingImpl(const std::string& parameter)
	: feature_names(),
	  attribute_names(),
	  bidirectional(false),
	  monotonicity(false),
	  lattice(0),
	  attr_phrase_span_first("phrase-span-first"),
	  attr_phrase_span_last("phrase-span-last")
      {
	typedef cicada::Parameter parameter_type;
      
	const parameter_type param(parameter);

	if (param.name() != "lexicalized-reordering"
	    && param.name() != "lexicalized-reorder"
	    && param.name() != "lexical-reordering"
	    && param.name() != "lexical-reorder")
	  throw std::runtime_error("is this really lexicalized reordering feature function? " + parameter);
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "bidirectional") == 0 || strcasecmp(piter->first.c_str(), "bi") == 0)
	    bidirectional = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "monotonicity") == 0 || strcasecmp(piter->first.c_str(), "mono") == 0)
	    monotonicity = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "feature") == 0)
	    attribute_names.push_back(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for lexicalized reordering: " << piter->first << "=" << piter->second << std::endl;
	}
	
	if (bidirectional) {
	  if (monotonicity) {
	    feature_names.reserve(4);
	    feature_names.push_back("lexicalized-reordering:forward-monotone");
	    feature_names.push_back("lexicalized-reordering:forward-others");
	    feature_names.push_back("lexicalized-reordering:backward-monotone");
	    feature_names.push_back("lexicalized-reordering:backward-others");
	  } else {
	    feature_names.reserve(6);
	    feature_names.push_back("lexicalized-reordering:forward-monotone");
	    feature_names.push_back("lexicalized-reordering:forward-swap");
	    feature_names.push_back("lexicalized-reordering:forward-discontinuous");
	    feature_names.push_back("lexicalized-reordering:backward-monotone");
	    feature_names.push_back("lexicalized-reordering:backward-swap");
	    feature_names.push_back("lexicalized-reordering:backward-discontinuous");
	  }
	} else {
	  if (monotonicity) {
	    feature_names.reserve(2);
	    feature_names.push_back("lexicalized-reordering:forward-monotone");
	    feature_names.push_back("lexicalized-reordering:forward-others");
	  } else {
	    feature_names.reserve(3);
	    feature_names.push_back("lexicalized-reordering:forward-monotone");
	    feature_names.push_back("lexicalized-reordering:forward-swap");
	    feature_names.push_back("lexicalized-reordering:forward-discontinuous");
	  }
	}
	
	attribute_names.resize(utils::bithack::max(feature_names.size(), attribute_names.size()));
	
	if (feature_names.size() != attribute_names.size())
	  throw std::runtime_error("attribute to feature mapping size do not match...");
	
	// assign default...
	for (size_t i = 0; i != attribute_names.size(); ++ i)
	  if (attribute_names[i].empty())
	    attribute_names[i] = "rule-table-" + boost::lexical_cast<std::string>(i);
      }

      void clear() { cache_states.clear(); }
      
      struct __feature_map : public boost::static_visitor<double>
      {
	double operator()(const double& x) const { return x; }
	template <typename Tp>
	double operator()(const Tp& x) const { return 0.0; }
      };
      
      void assign_feature(const feature_type& feature,
			  const double& score,
			  feature_set_type& features) const
      {
	if (score != 0.0)
	  features[feature] += score;
      }
      
      void reordering_score_next(const int prev_first, const int prev_last,
				 const int next_first, const int next_last,
				 const feature_cache_type& cache,
				 feature_set_type& features) const
      {
	if (bidirectional) {
	  if (monotonicity) {
	    if (prev_last == next_first)
	      assign_feature(feature_names[2], cache[2], features);
	    else
	      assign_feature(feature_names[3], cache[3], features);
	  } else {
	    if (prev_last == next_first)
	      assign_feature(feature_names[3], cache[3], features);
	    else if (next_last == prev_first)
	      assign_feature(feature_names[4], cache[4], features);
	    else
	      assign_feature(feature_names[5], cache[5], features);
	  }
	}
      }
      
      void reordering_score(const int prev_first, const int prev_last,
			    const int next_first, const int next_last,
			    const feature_cache_type& cache,
			    feature_set_type& features) const
      {
	if (monotonicity) {
	  if (prev_last == next_first)
	    assign_feature(feature_names[0], cache[0], features);
	  else
	    assign_feature(feature_names[1], cache[1], features);
	} else {
	  if (prev_last == next_first)
	    assign_feature(feature_names[0], cache[0], features);
	  else if (next_last == prev_first)
	    assign_feature(feature_names[1], cache[1], features);
	  else
	    assign_feature(feature_names[2], cache[2], features);
	}
      }
    
      void reordering_score_adjust(const int prev_first, const int prev_last,
				   const int next_first, const int next_last,
				   const feature_cache_type& cache,
				   feature_set_type& features) const
      {
	if (monotonicity) {
	  if (prev_last == next_first)
	    assign_feature(feature_names[0], - cache[0], features);
	  else
	    assign_feature(feature_names[1], - cache[1], features);
	} else {
	  if (prev_last == next_first)
	    assign_feature(feature_names[0], - cache[0], features);
	  else if (next_last == prev_first)
	    assign_feature(feature_names[1], - cache[1], features);
	  else
	    assign_feature(feature_names[2], - cache[2], features);
	}
      }
      
      struct __phrase_span : public boost::static_visitor<int>
      {
	int operator()(const attribute_set_type::int_type& x) const { return x; }
	template <typename Tp>
	int operator()(const Tp& x) const { throw std::runtime_error("no phrasal span with integer?"); }
      };
      
      int phrase_span(const attribute_set_type& attrs, const attribute_type& attr) const
      {
	attribute_set_type::const_iterator iter = attrs.find(attr);
	if (iter == attrs.end())
	  throw std::runtime_error("no phrasal span attribute?");
	
	return boost::apply_visitor(__phrase_span(), iter->second);
      }

      
      void lexicalized_reordering_score(state_ptr_type& state,
					const state_ptr_set_type& states,
					const edge_type& edge,
					feature_set_type& features) const
      {
	int* span = reinterpret_cast<int*>(state);
	state_type* node = reinterpret_cast<state_type*>(span + 2);
	
	if (states.empty()) {
	  // How do we capture initial phrase....???
	  const int span_first = phrase_span(edge.attributes, attr_phrase_span_first);
	  const int span_last  = phrase_span(edge.attributes, attr_phrase_span_last);

	  span[0] = span_first;
	  span[1] = span_last;
	  
	  feature_cache_type cache(feature_names.size(), 0.0);
	  for (size_t i = 0; i != feature_names.size(); ++ i) {
	    attribute_set_type::const_iterator aiter = edge.attributes.find(attribute_names[i]);
	    if (aiter != edge.attributes.end())
	      cache[i] = boost::apply_visitor(__feature_map(), aiter->second);
	  }
	  
	  feature_cache_states_type::iterator siter = const_cast<feature_cache_states_type&>(cache_states).insert(cache).first;
	  *node = siter - cache_states.begin();
	  
	  reordering_score(0, 0, span_first, span_last, cache, features);
	} else if (states.size() == 1) {
	  // it is only for the goal state...
	  const int*        span_antecedent = reinterpret_cast<const int*>(states[0]);
	  const state_type* node_antecedent = reinterpret_cast<const state_type*>(span_antecedent + 2);

	  span[0] = span_antecedent[0];
	  span[1] = span_antecedent[1];
	  *node   = *node_antecedent;
	} else if (states.size() == 2) {
	  const int*        span_antecedent = reinterpret_cast<const int*>(states[0]);
	  const state_type* node_antecedent = reinterpret_cast<const state_type*>(span_antecedent + 2);
	  
	  const int*        span_phrase     = reinterpret_cast<const int*>(states[1]);
	  const state_type* node_phrase     = reinterpret_cast<const state_type*>(span_phrase + 2);
	  
	  span[0] = span_phrase[0];
	  span[1] = span_phrase[1];
	  *node   = *node_phrase;
	  
	  // make adjustment, since this span-phrase is not a initial phrase..
	  reordering_score_adjust(0, 0, span_phrase[0], span_phrase[1], cache_states[*node_phrase], features);
	  reordering_score(span_antecedent[0], span_antecedent[1], span_phrase[0], span_phrase[1], cache_states[*node_phrase], features);
	  reordering_score_next(span_antecedent[0], span_antecedent[1], span_phrase[0], span_phrase[1], cache_states[*node_antecedent], features);
	} else
	  throw std::runtime_error("we do not support non-phrasal composed hypergraph");
      }
      
      void lexicalized_reordering_final_score(const state_ptr_type& state,
						feature_set_type& features) const
      {
	const int*        span = reinterpret_cast<const int*>(state);
	const state_type* node = reinterpret_cast<const state_type*>(span + 2);
	
	reordering_score_next(span[0], span[1], lattice->size(), lattice->size(), cache_states[*node], features);
      }
      
      void assign(const lattice_type& __lattice)
      {
	lattice = &__lattice;
      }
      
      feature_cache_states_type cache_states;
      
      feature_list_type   feature_names;
      attribute_list_type attribute_names;
      
      bool bidirectional;
      bool monotonicity;

      const lattice_type* lattice;

      attribute_type attr_phrase_span_first;
      attribute_type attr_phrase_span_last;
    };

    
    LexicalizedReordering::LexicalizedReordering(const std::string& parameter)
      : pimpl(0)
    {
      pimpl = new impl_type(parameter);
      
      // distotion context: span = [first, last)
      base_type::__state_size = sizeof(int) * 2 + sizeof(impl_type::state_type);
      base_type::__feature_name = std::string("lexicalized-reordering");
    }
    
    LexicalizedReordering::~LexicalizedReordering() { std::auto_ptr<impl_type> tmp(pimpl); }

    
    LexicalizedReordering::LexicalizedReordering(const LexicalizedReordering& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    { }
    
    
    LexicalizedReordering& LexicalizedReordering::operator=(const LexicalizedReordering& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void LexicalizedReordering::apply(state_ptr_type& state,
				      const state_ptr_set_type& states,
				      const edge_type& edge,
				      feature_set_type& features,
				      feature_set_type& estimates,
				      const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      pimpl->lexicalized_reordering_score(state, states, edge, features);
      
      if (final)
	pimpl->lexicalized_reordering_final_score(state, features);
    }
    
    void LexicalizedReordering::apply_coarse(state_ptr_type& state,
					     const state_ptr_set_type& states,
					     const edge_type& edge,
					     feature_set_type& features,
					     feature_set_type& estimates,
					     const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void LexicalizedReordering::apply_predict(state_ptr_type& state,
					      const state_ptr_set_type& states,
					      const edge_type& edge,
					      feature_set_type& features,
					      feature_set_type& estimates,
					      const bool final) const
    {}
    
    void LexicalizedReordering::apply_scan(state_ptr_type& state,
					   const state_ptr_set_type& states,
					   const edge_type& edge,
					   const int dot,
					   feature_set_type& features,
					   feature_set_type& estimates,
					   const bool final) const
    {}
    void LexicalizedReordering::apply_complete(state_ptr_type& state,
					       const state_ptr_set_type& states,
					       const edge_type& edge,
					       feature_set_type& features,
					       feature_set_type& estimates,
					       const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    void LexicalizedReordering::assign(const size_type& id,
				       const hypergraph_type& hypergraph,
				       const lattice_type& lattice,
				       const span_set_type& spans,
				       const sentence_set_type& targets,
				       const ngram_count_set_type& ngram_counts)
    {
      pimpl->assign(lattice);
    }

    void LexicalizedReordering::initialize()
    {
      pimpl->clear();
    }
  };
};
