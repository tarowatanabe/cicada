
#include <utility>
#include <memory>

#include "cicada/feature/boundary.hpp"
#include "cicada/parameter.hpp"

#include "utils/indexed_set.hpp"
#include "utils/compact_trie.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    
    class BoundaryImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef rule_type::symbol_set_type phrase_type;
      
      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;
      
      typedef utils::compact_trie<symbol_type, std::string, boost::hash<symbol_type>, std::equal_to<symbol_type>,
				  std::allocator<std::pair<const symbol_type, std::string> > > tree_map_type;

      typedef tree_map_type::id_type id_type;

      typedef std::pair<symbol_type, symbol_type> symbol_pair_type;
      typedef std::vector<symbol_pair_type, std::allocator<symbol_pair_type> > symbol_pair_set_type;

      typedef std::vector<int, std::allocator<int> > position_map_type;
      
      
      BoundaryImpl() : forced_feature(false) {}
      BoundaryImpl(const BoundaryImpl& x) : forced_feature(x.forced_feature) {}
      BoundaryImpl& operator=(const BoundaryImpl& x)
      {
	forced_feature = x.forced_feature;
	return *this;
      }

      virtual ~BoundaryImpl() {}
      
      void boundary_score(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features) const
      {
	if (states.empty()) {
	  symbol_type source_prefix = vocab_type::EPSILON;
	  symbol_type source_suffix = vocab_type::EPSILON;
	  symbol_type target_prefix = vocab_type::EPSILON;
	  symbol_type target_suffix = vocab_type::EPSILON;

	  compute_bound(edge.rule->source.begin(), edge.rule->source.end(), source_prefix, source_suffix);
	  compute_bound(edge.rule->target.begin(), edge.rule->target.end(), target_prefix, target_suffix);
	  
	  symbol_type* context = reinterpret_cast<symbol_type*>(state);
	  context[0] = source_prefix;
	  context[1] = source_suffix;
	  context[2] = target_prefix;
	  context[3] = target_suffix;
	} else {
	  phrase_span_set_type& source_spans = const_cast<phrase_span_set_type&>(source_spans_impl);
	  phrase_span_set_type& target_spans = const_cast<phrase_span_set_type&>(target_spans_impl);
	  
	  source_spans.clear();
	  target_spans.clear();
	  
	  edge.rule->source.terminals(std::back_inserter(source_spans));
	  edge.rule->target.terminals(std::back_inserter(target_spans));
	  
	  if (source_spans.size() != states.size() + 1)
	    throw std::runtime_error("# of states does not match...");
	  if (target_spans.size() != states.size() + 1)
	    throw std::runtime_error("# of states does not match...");
	  
	  // first, compute target side ranges and position mapping
	  
	  symbol_pair_set_type& symbol_pairs = const_cast<symbol_pair_set_type&>(symbol_pairs_impl);
	  position_map_type&    positions    = const_cast<position_map_type&>(positions_impl);
	  
	  symbol_pairs.clear();
	  positions.resize(states.size());

	  {
	    symbol_type prefix;
	    symbol_type suffix;
	    compute_bound(target_spans.front().first, target_spans.front().second, prefix, suffix);
	    
	    symbol_pairs.push_back(std::make_pair(prefix, suffix));
	    
	    phrase_span_set_type::const_iterator titer_begin = target_spans.begin();
	    phrase_span_set_type::const_iterator titer_end = target_spans.end();
	    for (phrase_span_set_type::const_iterator titer = titer_begin + 1; titer != titer_end; ++ titer) {
	      const phrase_span_type& span = *titer;
	      
	      // incase, we are working with non-synchronous parsing!
	      int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	      if (antecedent_index < 0)
		antecedent_index = titer - (titer_begin + 1);
	      
	      symbol_type prefix;
	      symbol_type suffix;
	      compute_bound(span.first, span.second, prefix, suffix);
	      
	      symbol_pairs.push_back(std::make_pair(prefix, suffix));
	      
	      positions[antecedent_index] = titer - (titer_begin + 1);
	    }
	  }
	  
	  // we keep source-prefix and source-suffix...
	  symbol_type prefix;
	  symbol_type suffix;
	  
	  compute_bound(source_spans.front().first, source_spans.front().second, prefix, suffix);

	  phrase_span_set_type::const_iterator siter_begin = source_spans.begin();
	  phrase_span_set_type::const_iterator siter_end = source_spans.end();
	  for (phrase_span_set_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    // incase, we are working with non-synchronous parsing!
	    int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      antecedent_index = siter - (siter_begin + 1);
	    
	    const symbol_type* antecedent_context = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	    const symbol_type& antecedent_source_prefix = antecedent_context[0];
	    const symbol_type& antecedent_source_suffix = antecedent_context[1];
	    const symbol_type& antecedent_target_prefix = antecedent_context[2];
	    const symbol_type& antecedent_target_suffix = antecedent_context[3];
	    
	    const symbol_type target_suffix      = boundary_suffix(positions[antecedent_index], states, symbol_pairs);
	    const symbol_type target_prefix_next = boundary_prefix(positions[antecedent_index] + 1, states, symbol_pairs);
	    
	    symbol_type prefix_next;
	    symbol_type suffix_next;
	    compute_bound(span.first, span.second, prefix_next, suffix_next);
	    
	    if (! suffix.empty() && ! target_suffix.empty())
	      apply_feature(features, suffix, antecedent_source_prefix, target_suffix, antecedent_target_prefix);
	    
	    if (prefix.empty())
	      prefix = antecedent_source_prefix;
	    
	    if (! prefix_next.empty()) {
	      if (! target_prefix_next.empty())
		apply_feature(features, antecedent_source_suffix, prefix_next, antecedent_target_suffix, target_prefix_next);
	      
	      suffix = suffix_next;
	    } else
	      suffix = antecedent_source_suffix;
	  }
	  
	  symbol_type* context = reinterpret_cast<symbol_type*>(state);
	  context[0] = (prefix.empty() ? vocab_type::EPSILON : prefix);
	  context[1] = (suffix.empty() ? vocab_type::EPSILON : suffix);
	  context[2] = boundary_prefix(0, states, symbol_pairs);
	  context[3] = boundary_suffix(symbol_pairs.size() - 1, states, symbol_pairs);
	}
      }
      
      void boundary_final_score(const state_ptr_type& state,
				feature_set_type& features) const
      {
	const symbol_type* context = reinterpret_cast<const symbol_type*>(state);
	const symbol_type& source_prefix = context[0];
	const symbol_type& source_suffix = context[1];
	const symbol_type& target_prefix = context[2];
	const symbol_type& target_suffix = context[3];
	
	apply_feature(features, vocab_type::BOS, source_prefix,   vocab_type::BOS, target_prefix);
	apply_feature(features, source_suffix,   vocab_type::EOS, target_suffix,   vocab_type::EOS);
      }
      
      void apply_feature(feature_set_type& features,
			 const std::string& source_prev, const std::string& source_next,
			 const std::string& target_prev, const std::string& target_next) const
      {
	const std::string name = "boundary:" + source_prev + '|' + source_next + '|' + target_prev + '|' + target_next;
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
      }
      
      const symbol_type& boundary_prefix(const int index,
					 const state_ptr_set_type& states,
					 const symbol_pair_set_type& symbol_pairs) const
      {
	const symbol_type& prefix = symbol_pairs[index].first;

	if (index + 1 == symbol_pairs.size() || ! prefix.empty())
	  return prefix;
	else
	  return reinterpret_cast<const symbol_type*>(states[index])[2];
      }

      const symbol_type& boundary_suffix(const int index,
					 const state_ptr_set_type& states,
					 const symbol_pair_set_type& symbol_pairs) const
      {
	const symbol_type& suffix = symbol_pairs[index].second;
	
	if (index == 0 || ! suffix.empty())
	  return suffix;
	else
	  return reinterpret_cast<const symbol_type*>(states[index - 1])[3];
      }
      
      template <typename Iterator>
      void compute_bound(Iterator first, Iterator last, symbol_type& prefix, symbol_type& suffix) const
      {
	for (Iterator iter = first; iter != last; ++ iter)
	  if (*iter != vocab_type::EPSILON) {
	    prefix = *iter;
	    break;
	  }
	
	for (Iterator iter = last; iter != first; -- iter)
	  if (*(iter - 1) != vocab_type::EPSILON) {
	    suffix = *(iter - 1);
	    break;
	  }
      }
      
      phrase_span_set_type source_spans_impl;
      phrase_span_set_type target_spans_impl;

      symbol_pair_set_type symbol_pairs_impl;
      position_map_type    positions_impl;

      bool forced_feature;
    };
    
    
    Boundary::Boundary(const std::string& parameter)
      : pimpl(new impl_type())
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (param.name() != "boundary")
	throw std::runtime_error("is this really boundary feature function? " + parameter);
      
      base_type::__state_size = sizeof(symbol_type) * 4;
      base_type::__feature_name = "boundary";
      base_type::__sparse_feature = true;
    }
    
    Boundary::~Boundary() { std::auto_ptr<impl_type> tmp(pimpl); }

    Boundary::Boundary(const Boundary& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    Boundary& Boundary::operator=(const Boundary& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    template <typename FeaturePrefix, typename Feature>
    inline
    bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
    {
      return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
    }
    
    void Boundary::operator()(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      feature_set_type& features,
			      feature_set_type& estimates) const
    {
      const std::string& __feature_prefix = base_type::feature_name();
      for (feature_set_type::iterator fiter = features.begin(); fiter != features.end(); /**/)
	if (equal_prefix(__feature_prefix, fiter->first))
	  features.erase(fiter ++);
	else
	  ++ fiter;
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->boundary_score(state, states, edge, features);
    }
    
    void Boundary::operator()(const state_ptr_type& state,
			      feature_set_type& features,
			      feature_set_type& estimates) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->boundary_final_score(state, features);
    }

    void Boundary::initialize()
    {
      
    }
  };
};
