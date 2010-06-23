
#include <utility>
#include <memory>

#include "cicada/feature/neighbours.hpp"
#include "cicada/parameter.hpp"

#include "utils/indexed_set.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    
    class NeighboursImpl
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

      struct state_type
      {
	symbol_type node;
	symbol_type prefix;
	symbol_type suffix;
	int         span;

	state_type()
	  : node(), prefix(), suffix(), span() {}

	state_type(const symbol_type& __node,
		   const symbol_type& __prefix,
		   const symbol_type& __suffix,
		   const int&         __span)
	  : node(__node), prefix(__prefix), suffix(__suffix), span(__span) {}

	friend
	bool operator==(const state_type& x, const state_type& y)
	{
	  return x.node == y.node && x.prefix == y.prefix && x.suffix == y.suffix && x.span == y.span;
	}
      };
      
      struct state_hash_type : public utils::hashmurmur<size_t>
      {
	typedef utils::hashmurmur<size_t> hasher_type;
	
	size_t operator()(const state_type& state) const
	{
	  return hasher_type::operator()(state, 0);
	}
      };
      
      typedef utils::indexed_set<state_type, state_hash_type, std::equal_to<state_type>, std::allocator<state_type> > state_set_type;
      
      typedef uint32_t id_type;
      
      virtual ~NeighboursImpl() {}
      
      virtual void neighbours_score(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features) const = 0;
      virtual void neighbours_final_score(const state_ptr_type& state,
					  feature_set_type& features) const = 0;
      
      void clear()
      {
	
	
      }
      
      state_set_type states_id;
      
      phrase_span_set_type phrase_spans_impl;
    };
    
    template <typename Extract>
    class __NeighboursImpl : public NeighboursImpl, public Extract
    {
      
      virtual void neighbours_score(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features) const
      {
	// this feature function is complicated in that we know nothing about the source-side...

	const rule_type::symbol_set_type& phrase = extract_phrase(edge);
	
	if (states.empty()) {
	  // we do not add feature here, since we know nothing abount surrounding context...
	  symbol_type* context_node   = reinterpret_cast<symbol_type*>(state);
	  symbol_type* context_prefix = context_node + 1;
	  symbol_type* context_suffix = context_prefix + 1;
	  int* context_size = reinterpret_cast<int*>(context_suffix + 1);
	  
	  symbol_type prefix = vocab_type::EPSILON;
	  symbol_type suffix = vocab_type::EPSILON;
	  int size = 0;
	  
	  rule_type::symbol_set_type::const_iterator piter_end = phrase.end();
	  for (rule_type::symbol_set_type::const_iterator piter = phrase.begin(); piter != piter_end; ++ piter) 
	    if (*piter != vocab_type::EPSILON) {
	      if (size == 0)
		prefix = *piter;
	      suffix = *piter;
	      ++ size;
	    }
	  
	  *context_node   = edge.rule->lhs;
	  *context_prefix = prefix;
	  *context_suffix = suffix;
	  *context_size   = size;

	} else {
	  phrase_span_set_type& phrase_spans = const_cast<phrase_span_set_type&>(phrase_spans_impl);
	  
	  phrase_spans.clear();
	  phrase.terminals(std::back_inserter(phrase_spans));

	  if (phrase_spans.size() != states.size() + 1)
	    throw std::runtime_error("# of states does not match...");

	  symbol_type prefix = vocab_type::EPSILON;
	  symbol_type suffix = vocab_type::EPSILON;
	  int size = count_span(phrase_spans.front().first, phrase_spans.front().second, prefix, suffix);
	  
	  phrase_span_set_type::const_iterator siter_begin = phrase_spans.begin();
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    // incase, we are working with non-synchronous parsing!
	    if (antecedent_index < 0)
	      antecedent_index = siter - (siter_begin + 1);
	    
	    const symbol_type* context_node   = reinterpret_cast<const symbol_type*>(states[antecedent_index]);
	    const symbol_type* context_prefix = context_node + 1;
	    const symbol_type* context_suffix = context_prefix + 1;
	    const int* context_size = reinterpret_cast<const int*>(context_suffix + 1);
	    
	    symbol_type prefix_next = vocab_type::EPSILON;
	    symbol_type suffix_next = vocab_type::EPSILON;
	    int size_next = count_span(span.first, span.second, prefix_next, suffix_next);

	    if (suffix != vocab_type::EPSILON) {
	    
	      if (prefix_next != vocab_type::EPSILON) {
		const std::string& node = static_cast<const std::string&>(*context_node);
		const std::string& prev = static_cast<const std::string&>(suffix);
		const std::string& next = static_cast<const std::string&>(prefix_next);
		
		// perform scoring for this antecedent
		features[feature_prefix() + node + '|' + prev + '|' + next + '|' + boost::lexical_cast<std::string>(*context_size)] += 1.0;
		
		suffix = suffix_next;
		size += *context_size + size_next;
	      } else if (siter + 1 != siter_end) {
		// we have prefix and suffix from next antecedent node...
		
	      } else {
		// we have prefix, but no suffix...
		
	      }
	    } else {
	      if (siter - 1 != siter_begin) {
		// we have suffix and prefix from previous antecedent node
		
	      } else {
		// we have no
	      } 
	    }
	  }
	}
      }
      
      virtual void neighbours_final_score(const state_ptr_type& state,
					  feature_set_type& features) const
      {
	const symbol_type* context_node   = reinterpret_cast<const symbol_type*>(state);
	const symbol_type* context_prefix = context_node + 1;
	const symbol_type* context_suffix = context_prefix + 1;
	const int* context_size = reinterpret_cast<const int*>(context_suffix + 1);
	
	const std::string& bos  = static_cast<const std::string&>(vocab_type::BOS);
	const std::string& eos  = static_cast<const std::string&>(vocab_type::EOS);
	
	features[feature_prefix() + static_cast<const std::string&>(*context_node) + '|' + bos + '|' + eos + '|' + boost::lexical_cast<std::string>(*context_size)] += 1.0;
      }
	  
      template <typename Iterator>
      int count_span(Iterator first, Iterator last, symbol_type& prefix, symbol_type& suffix) const
      {
	int count = 0;
	for (/**/; first != last; ++ first) 
	  if (*first != vocab_type::EPSILON) {
	    if (count == 0)
	      prefix = *first;
	    suffix = *first;
	    ++ count;
	  }
	return count;
      }

      const std::string& feature_prefix() const
      {
	return Extract::feature_prefix;
      }

      template <typename Edge>
      const rule_type::symbol_set_type& extract_phrase(const Edge& x) const
      {
	static const rule_type::symbol_set_type __tmptmp;

	return Extract::operator()(x, __tmptmp);
      }
    };
    

    struct __neighbours_extract_source
    {
      __neighbours_extract_source()
	: feature_prefix("neighbours-source:") {}
      
      template <typename Edge, typename Phrase>
      const Phrase& operator()(const Edge& x, const Phrase& phrase) const
      {
	return x.rule->source;
      }
      
      const std::string feature_prefix;
    };

    struct __neighbours_extract_target
    {
      __neighbours_extract_target()
	: feature_prefix("neighbours-target:") {}
      
      template <typename Edge, typename Phrase>
      const Phrase& operator()(const Edge& x, const Phrase& phrase) const
      {
	return x.rule->target;
      }
      
      const std::string feature_prefix;
    };
    
    inline
    bool true_false(const std::string& token)
    {
      if (strcasecmp(token.c_str(), "true") == 0)
	return true;
      if (strcasecmp(token.c_str(), "yes") == 0)
	return true;
      if (atoi(token.c_str()) > 0)
	return true;
      return false;
    }

    
    Neighbours::Neighbours(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "neighbours")
	throw std::runtime_error("is this really neighbours feature function? " + parameter);

      bool source = false;
      bool target = false;
      if (param.find("source") != param.end())
	source = true_false(param.find("source")->second);

      if (param.find("target") != param.end())
	target = true_false(param.find("target")->second);
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      if (! source && ! target)
	throw std::runtime_error("what side are you going to use?");
      
      std::auto_ptr<impl_type> neighbours_impl(source
					       ? dynamic_cast<impl_type*>(new __NeighboursImpl<__neighbours_extract_source>())
					       : dynamic_cast<impl_type*>(new __NeighboursImpl<__neighbours_extract_target>()));

      
      // non-terminal + two neighbouring symbols + span-size
      base_type::__state_size = sizeof(symbol_type) * 3 + sizeof(int);
      
      pimpl = neighbours_impl.release();
    }
    
    Neighbours::~Neighbours() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    void Neighbours::operator()(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates) const
    {
      pimpl->neighbours_score(state, states, edge, features);
    }
    
    void Neighbours::operator()(const state_ptr_type& state,
				feature_set_type& features) const
    {
      pimpl->neighbours_final_score(state, features);
    }

  };
};
