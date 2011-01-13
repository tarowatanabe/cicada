//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>
#include <algorithm>

#include "cicada/feature/neighbours.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"
#include "cicada/cluster_stemmer.hpp"

#include "utils/indexed_set.hpp"
#include "utils/lexical_cast.hpp"

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

      typedef cicada::Cluster  cluster_type;
      typedef cicada::Stemmer  stemmer_type;

      typedef cicada::ClusterStemmer normalizer_type;
      typedef std::vector<normalizer_type, std::allocator<normalizer_type> > normalizer_set_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef rule_type::symbol_set_type phrase_type;
      
      typedef std::pair<phrase_type::const_iterator, phrase_type::const_iterator> phrase_span_type;
      typedef std::vector<phrase_span_type, std::allocator<phrase_span_type> >  phrase_span_set_type;

      typedef uint32_t id_type;

      struct state_type
      {
	id_type     parent;
	symbol_type node;
	symbol_type prefix;
	symbol_type suffix;
	int         span;
	
	state_type()
	  : parent(id_type(-1)), node(), prefix(), suffix(), span() {}
	
	state_type(const id_type&     __parent,
		   const symbol_type& __node,
		   const symbol_type& __prefix,
		   const symbol_type& __suffix,
		   const int&         __span)
	  : parent(__parent), node(__node), prefix(__prefix), suffix(__suffix), span(__span) {}

	friend
	bool operator==(const state_type& x, const state_type& y)
	{
	  return x.parent == y.parent && x.node == y.node && x.prefix == y.prefix && x.suffix == y.suffix && x.span == y.span;
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
      
      typedef utils::indexed_set<state_type, state_hash_type, std::equal_to<state_type>, std::allocator<state_type> > state_map_type;
      
      NeighboursImpl()
	: sentence(0),
	  forced_feature(false),
	  alignment_mode(false),
	  attr_target_position("target-position") {}
      
      void clear()
      {
	state_map.clear();
      }

      normalizer_set_type normalizers;
      
      state_map_type state_map;
      
      phrase_span_set_type phrase_spans_impl;
      
      feature_type feature_name_prefix;
      const sentence_type* sentence;
      
      bool forced_feature;
      bool alignment_mode;

      attribute_type attr_target_position;

      struct __attribute_integer : public boost::static_visitor<cicada::AttributeVector::int_type>
      {
	typedef cicada::AttributeVector attribute_set_type;
	
	attribute_set_type::int_type operator()(const attribute_set_type::int_type& x) const { return x; }
	attribute_set_type::int_type operator()(const attribute_set_type::float_type& x) const { return -2; }
	attribute_set_type::int_type operator()(const attribute_set_type::string_type& x) const { return -2; }
      };
      
      
      void neighbours_score(state_ptr_type& state,
			    const state_ptr_set_type& states,
			    const edge_type& edge,
			    feature_set_type& features) const
      {
	// this feature function is complicated in that we know nothing about the source-side...

	const rule_type::symbol_set_type& phrase = edge.rule->rhs;
	
	if (states.empty()) {
	  // we do not add feature here, since we know nothing abount surrounding context...
	  
	  if (alignment_mode) {
	    symbol_type prefix = vocab_type::EPSILON;
	    symbol_type suffix = vocab_type::EPSILON;
	    int size = 0;
	    
	    attribute_set_type::const_iterator titer = edge.attributes.find(attr_target_position);
	    if (titer == edge.attributes.end())
	      throw std::runtime_error("we do not support non alignment forest");
	    
	    const int target_pos = boost::apply_visitor(__attribute_integer(), titer->second);
	    
	    if (sentence && target_pos >= 0) {
	      prefix = sentence->operator[](target_pos);
	      suffix = sentence->operator[](target_pos);
	      size = 1;
	    }
	    
	    state_map_type::iterator iter = const_cast<state_map_type&>(state_map).insert(state_type(id_type(-1), edge.rule->lhs, prefix, suffix, size)).first;
	    
	    *reinterpret_cast<id_type*>(state) = iter - state_map.begin();
	  } else {
	    symbol_type prefix = vocab_type::EPSILON;
	    symbol_type suffix = vocab_type::EPSILON;
	    int size = count_span(phrase.begin(), phrase.end(), prefix, suffix);
	    
	    state_map_type::iterator iter = const_cast<state_map_type&>(state_map).insert(state_type(id_type(-1), edge.rule->lhs, prefix, suffix, size)).first;
	    
	    *reinterpret_cast<id_type*>(state) = iter - state_map.begin();
	  }
	} else {
	  phrase_span_set_type& phrase_spans = const_cast<phrase_span_set_type&>(phrase_spans_impl);
	  
	  phrase_spans.clear();
	  phrase.terminals(std::back_inserter(phrase_spans));
	  
	  if (phrase_spans.size() != states.size() + 1)
	    throw std::runtime_error("# of states does not match...");
	  
	  symbol_type prefix = vocab_type::EMPTY;
	  symbol_type suffix = vocab_type::EMPTY;
	  int size = count_span(phrase_spans.front().first, phrase_spans.front().second, prefix, suffix);

	  id_type state_id(id_type(-1));
	  
	  phrase_span_set_type::const_iterator siter_begin = phrase_spans.begin();
	  phrase_span_set_type::const_iterator siter_end = phrase_spans.end();
	  for (phrase_span_set_type::const_iterator siter = siter_begin + 1; siter != siter_end; ++ siter) {
	    const phrase_span_type& span = *siter;
	    
	    // incase, we are working with non-synchronous parsing!
	    int antecedent_index = (span.first - 1)->non_terminal_index() - 1;
	    if (antecedent_index < 0)
	      antecedent_index = siter - (siter_begin + 1);
	    
	    const id_type antecedent_id = *reinterpret_cast<const id_type*>(states[antecedent_index]);

	    const symbol_type prefix_antecedent = state_map[antecedent_id].prefix;
	    const symbol_type suffix_antecedent = state_map[antecedent_id].suffix;
	    
	    symbol_type prefix_next = vocab_type::EMPTY;
	    symbol_type suffix_next = vocab_type::EMPTY;
	    const int size_next = count_span(span.first, span.second, prefix_next, suffix_next);

	    size  += size_next + state_map[antecedent_id].span;
	    
	    if (prefix == vocab_type::EMPTY)
	      prefix = prefix_antecedent;
	    
	    if (prefix_next != vocab_type::EMPTY) {
	      state_id = apply_features(features, state_id, antecedent_id, suffix, prefix_next);
	      
	      suffix = suffix_next;
	    } else if (siter + 1 != siter_end) {
	      // we have no prefix for this span, but prefix from next antecedent node...
	      //thus, use next antecedent's prefix as our suffix
	      
	      int antecedent_index_next = (span.second)->non_terminal_index() - 1;
	      if (antecedent_index_next < 0)
		antecedent_index_next = siter + 1 - (siter_begin + 1);
	      
	      const id_type antecedent_id_next = *reinterpret_cast<const id_type*>(states[antecedent_index_next]);
	      
	      state_id = apply_features(features, state_id, antecedent_id, suffix, state_map[antecedent_id_next].prefix);
	      
	      suffix = suffix_antecedent;
	    } else {
	      // we have nothing as our suffix
	      state_id = apply_features(features, state_id, antecedent_id, suffix, vocab_type::EMPTY);
	      
	      suffix = suffix_antecedent;
	    } 
	  }
	  
	  // construct state is used as supplier for the next state...
	  state_map_type::iterator iter = const_cast<state_map_type&>(state_map).insert(state_type(state_id, edge.rule->lhs, prefix, suffix, size)).first;
	  
	  *reinterpret_cast<id_type*>(state) = iter - state_map.begin();
	}
      }

      void apply_feature(feature_set_type& features, const std::string& node, const symbol_type& prev, const symbol_type& next, const int span) const
      {
	const std::string name = feature_name(node, prev, next, span);
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;

	for (size_t i = 0; i != normalizers.size(); ++ i) {
	  const symbol_type prev_norm = normalizers[i](prev);
	  const symbol_type next_norm = normalizers[i](next);
	  
	  if (prev_norm != prev || next_norm != next) {
	    const std::string name = feature_name(node, prev_norm, next_norm, span);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}
      }

      id_type apply_features(feature_set_type& features, id_type id_curr, id_type id, const symbol_type& prefix, const symbol_type& suffix) const
      {
	typedef std::vector<state_type, std::allocator<state_type> > state_set_type;

	state_set_type states;

	if (prefix == vocab_type::EMPTY && suffix == vocab_type::EMPTY) {
	  while (id != id_type(-1)) {
	    const state_type& state = state_map[id];
	    
	    if (state.prefix == vocab_type::EMPTY && state.suffix == vocab_type::EMPTY)
	      states.push_back(state_type(id_curr, state.node, prefix, suffix, state.span));
	    else if (state.prefix == vocab_type::EMPTY)
	      states.push_back(state_type(id_curr, state.node, prefix, state.suffix, state.span));
	    else if (state.suffix == vocab_type::EMPTY)
	      states.push_back(state_type(id_curr, state.node, state.prefix, suffix, state.span));
	    else
	      states.push_back(state_type(id_curr, state.node, prefix, suffix, state.span));
	    
	    id = state.parent;
	  }
	  
	} else if (prefix == vocab_type::EMPTY) {
	  // we have at least suffix context...
	  
	  while (id != id_type(-1)) {
	    const state_type& state = state_map[id];
	    
	    if (state.prefix == vocab_type::EMPTY && state.suffix == vocab_type::EMPTY)
	      states.push_back(state_type(id_curr, state.node, prefix, suffix, state.span));
	    else if (state.prefix == vocab_type::EMPTY)
	      states.push_back(state_type(id_curr, state.node, prefix, state.suffix, state.span));
	    else if (state.suffix == vocab_type::EMPTY)
	      apply_feature(features, state.node, state.prefix, suffix, state.span);
	    else
	      states.push_back(state_type(id_curr, state.node, prefix, suffix, state.span));
	    
	    id = state.parent;
	  }
	} else if (suffix == vocab_type::EMPTY) {
	  // we have at least prefix context...
	  
	  while (id != id_type(-1)) {
	    const state_type& state = state_map[id];
	    
	    if (state.prefix == vocab_type::EMPTY && state.suffix == vocab_type::EMPTY)
	      states.push_back(state_type(id_curr, state.node, prefix, suffix, state.span));
	    else if (state.suffix == vocab_type::EMPTY)
	      states.push_back(state_type(id_curr, state.node, state.prefix, suffix, state.span));
	    else if (state.prefix == vocab_type::EMPTY)
	      apply_feature(features, state.node, prefix, state.suffix, state.span);
	    else
	      states.push_back(state_type(id_curr, state.node, prefix, suffix, state.span));
	    
	    id = state.parent;
	  }
	} else {
	  while (id != id_type(-1)) {
	    const state_type& state = state_map[id];
	  
	    if (state.prefix == vocab_type::EMPTY && state.suffix == vocab_type::EMPTY)
	      apply_feature(features, state.node, prefix, suffix, state.span);
	    else if (state.prefix == vocab_type::EMPTY)
	      apply_feature(features, state.node, prefix, state.suffix, state.span);
	    else if (state.suffix == vocab_type::EMPTY)
	      apply_feature(features, state.node, state.prefix, suffix, state.span);
	    else
	      apply_feature(features, state.node, prefix, suffix, state.span);
	  
	    id = state.parent;
	  }
	}
	
	state_set_type::const_reverse_iterator siter_end = states.rend();
	for (state_set_type::const_reverse_iterator siter = states.rbegin(); siter != siter_end; ++ siter) {
	  state_map_type::iterator iter = const_cast<state_map_type&>(state_map).insert(state_type(id_curr, siter->node, siter->prefix, siter->suffix, siter->span)).first;
	  id_curr = iter - state_map.begin();
	}
	
	return id_curr;
      }
      
      const std::string feature_name(const std::string& node, const std::string& prev, const std::string& next, const int span) const
      {
	return (static_cast<const std::string&>(feature_name_prefix) + ":" + node + '|' + prev + '|' + next + '|' + boost::lexical_cast<std::string>(span));
      }
      
      void neighbours_final_score(const state_ptr_type& __state,
				  feature_set_type& features) const
      {

	apply_features(features, id_type(-1), *reinterpret_cast<const id_type*>(__state), vocab_type::BOS, vocab_type::EOS);
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
    };
    
    
    Neighbours::Neighbours(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "neighbours" && param.name() != "neighbors")
	throw std::runtime_error("is this really neighbours feature function? " + parameter);

      impl_type::normalizer_set_type normalizers;
      std::string name;
      bool alignment_mode = false;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "cluster") == 0) {
	  if (! boost::filesystem::exists(piter->second))
	    throw std::runtime_error("no cluster file: " + piter->second);
	  
	  normalizers.push_back(impl_type::normalizer_type(&cicada::Cluster::create(piter->second)));
	} else if (strcasecmp(piter->first.c_str(), "stemmer") == 0)
	  normalizers.push_back(impl_type::normalizer_type(&cicada::Stemmer::create(piter->second)));
	else if (strcasecmp(piter->first.c_str(), "name") == 0)
	  name = piter->second;
	else if (strcasecmp(piter->first.c_str(), "alignment") == 0)
	  alignment_mode = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for neighbours: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> neighbours_impl(new impl_type());
      
      neighbours_impl->normalizers.swap(normalizers);
      neighbours_impl->alignment_mode = alignment_mode;
      neighbours_impl->feature_name_prefix = (name.empty() ? std::string("neighbours") : name);
      
      // non-terminal + two neighbouring symbols + span-size
      base_type::__state_size = sizeof(impl_type::id_type);
      base_type::__feature_name = (name.empty() ? std::string("neighbours") : name);
      base_type::__sparse_feature = true;
      
      pimpl = neighbours_impl.release();
    }
    
    Neighbours::~Neighbours() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    Neighbours::Neighbours(const Neighbours& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    { }
    
    Neighbours& Neighbours::operator=(const Neighbours& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      return *this;
    }
    
    void Neighbours::apply(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   feature_set_type& features,
			   feature_set_type& estimates,
			   const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->neighbours_score(state, states, edge, features);

      if (final)
	pimpl->neighbours_final_score(state, features);
    }

    void Neighbours::apply_coarse(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {}
    
    void Neighbours::apply_predict(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {}
    
    void Neighbours::apply_scan(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				const int dot,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {}
    
    void Neighbours::apply_complete(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    void Neighbours::initialize()
    {
      pimpl->clear();
    }

    
    void Neighbours::assign(const size_type& id,
			    const hypergraph_type& hypergraph,
			    const lattice_type& lattice,
			    const span_set_type& spans,
			    const sentence_set_type& targets,
			    const ngram_count_set_type& ngram_counts)
    {
      pimpl->sentence = 0;
      if (! targets.empty())
	pimpl->sentence = &targets.front();
    }

  };
};
