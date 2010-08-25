
#include <utility>
#include <memory>
#include <algorithm>

#include "cicada/feature/neighbours.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"

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
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type feature_set_type;
      
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
	: cluster(0), stemmer_prefix(0), stemmer_suffix(0), stemmer_digits(0), forced_feature(false) {}
      
      virtual ~NeighboursImpl() {}
      
      virtual void neighbours_score(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features) const = 0;
      virtual void neighbours_final_score(const state_ptr_type& state,
					  feature_set_type& features) const = 0;
      
      void clear()
      {
	state_map.clear();
      }

      cluster_type* cluster;
      stemmer_type* stemmer_prefix;
      stemmer_type* stemmer_suffix;
      stemmer_type* stemmer_digits;
      
      state_map_type state_map;
      
      phrase_span_set_type phrase_spans_impl;

      bool forced_feature;
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
	  symbol_type prefix = vocab_type::EPSILON;
	  symbol_type suffix = vocab_type::EPSILON;
	  int size = count_span(phrase.begin(), phrase.end(), prefix, suffix);
	  
	  state_map_type::iterator iter = const_cast<state_map_type&>(state_map).insert(state_type(id_type(-1), edge.rule->lhs, prefix, suffix, size)).first;

	  *reinterpret_cast<id_type*>(state) = iter - state_map.begin();
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
	      const state_type& state_next = state_map[antecedent_id_next];
	      
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
	if (cluster) {
	  const symbol_type prev_cluster = cluster->operator[](prev);
	  const symbol_type next_cluster = cluster->operator[](next);
	  
	  if (prev_cluster != prev || next_cluster != next) {
	    const std::string name = feature_name(node, prev_cluster, next_cluster, span);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}

	if (stemmer_prefix) {
	  const symbol_type prev_stemmed = stemmer_prefix->operator[](prev);
	  const symbol_type next_stemmed = stemmer_prefix->operator[](next);
	  
	  if (prev_stemmed != prev || next_stemmed != next) {
	    const std::string name = feature_name(node, prev_stemmed, next_stemmed, span);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}

	if (stemmer_suffix) {
	  const symbol_type prev_stemmed = stemmer_suffix->operator[](prev);
	  const symbol_type next_stemmed = stemmer_suffix->operator[](next);
	  
	  if (prev_stemmed != prev || next_stemmed != next) {
	    const std::string name = feature_name(node, prev_stemmed, next_stemmed, span);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}

	if (stemmer_digits) {
	  const symbol_type prev_stemmed = stemmer_digits->operator[](prev);
	  const symbol_type next_stemmed = stemmer_digits->operator[](next);
	  
	  if (prev_stemmed != prev || next_stemmed != next) {
	    const std::string name = feature_name(node, prev_stemmed, next_stemmed, span);
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}
	
	const std::string name = feature_name(node, prev, next, span);
	if (forced_feature || feature_set_type::feature_type::exists(name))
	  features[name] += 1.0;
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
	return Extract::feature_prefix + node + '|' + prev + '|' + next + '|' + boost::lexical_cast<std::string>(span);
      }
      
      virtual void neighbours_final_score(const state_ptr_type& __state,
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
      
      std::string feature_prefix;
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
      
      std::string feature_prefix;
    };
    
    Neighbours::Neighbours(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "neighbours" && param.name() != "neighbors")
	throw std::runtime_error("is this really neighbours feature function? " + parameter);
      
      bool source = false;
      bool target = false;
      
      int stemmer_prefix_size = 0;
      int stemmer_suffix_size = 0;
      bool stemmer_digits = false;
      
      boost::filesystem::path cluster_path;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	  const std::string& yield = piter->second;
	  
	  if (strcasecmp(yield.c_str(), "source") == 0)
	    source = true;
	  else if (strcasecmp(yield.c_str(), "target") == 0)
	    target = true;
	  else
	    throw std::runtime_error("unknown parameter: " + parameter);
	  
	} else if (strcasecmp(piter->first.c_str(), "cluster") == 0)
	  cluster_path = piter->second;
	else if (strcasecmp(piter->first.c_str(), "prefix") == 0)
	  stemmer_prefix_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "suffix") == 0)
	  stemmer_suffix_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "digits") == 0)
	  stemmer_digits = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for neighbours: " << piter->first << "=" << piter->second << std::endl;
      }
      
      if (source && target)
	throw std::runtime_error("both source and target?");
      if (! source && ! target)
	throw std::runtime_error("what side are you going to use?");
      
      if (stemmer_prefix_size < 0)
	throw std::runtime_error("negative prefix size?");
      if (stemmer_suffix_size < 0)
	throw std::runtime_error("negative suffix size?");

      std::auto_ptr<impl_type> neighbours_impl(source
					       ? dynamic_cast<impl_type*>(new __NeighboursImpl<__neighbours_extract_source>())
					       : dynamic_cast<impl_type*>(new __NeighboursImpl<__neighbours_extract_target>()));

      if (! cluster_path.empty()) {
	if (! boost::filesystem::exists(cluster_path))
	  throw std::runtime_error("no cluster file: " + cluster_path.file_string());

	neighbours_impl->cluster = &cicada::Cluster::create(cluster_path);
      }

      if (stemmer_prefix_size > 0)
	neighbours_impl->stemmer_prefix = &cicada::Stemmer::create("prefix:size=" + boost::lexical_cast<std::string>(stemmer_prefix_size));
      
      if (stemmer_suffix_size > 0)
	neighbours_impl->stemmer_suffix = &cicada::Stemmer::create("suffix:size=" + boost::lexical_cast<std::string>(stemmer_suffix_size));

      if (stemmer_digits)
	neighbours_impl->stemmer_digits = &cicada::Stemmer::create("digits");

      
      // non-terminal + two neighbouring symbols + span-size
      base_type::__state_size = sizeof(impl_type::id_type);
      base_type::__feature_name = std::string("neighbours-") + (source ? "source" : "target");
      base_type::__sparse_feature = true;
      
      pimpl = neighbours_impl.release();
    }
    
    Neighbours::~Neighbours() { std::auto_ptr<impl_type> tmp(pimpl); }
    
    Neighbours::Neighbours(const Neighbours& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(0)
    {
      typedef __NeighboursImpl<__neighbours_extract_source> neighbours_source_type;
      typedef __NeighboursImpl<__neighbours_extract_target> neighbours_target_type;
      
      if (dynamic_cast<const neighbours_source_type*>(x.pimpl))
	pimpl = new neighbours_source_type();
      else
	pimpl = new neighbours_target_type();
    }
    
    Neighbours& Neighbours::operator=(const Neighbours& x)
    {
      typedef __NeighboursImpl<__neighbours_extract_source> neighbours_source_type;
      typedef __NeighboursImpl<__neighbours_extract_target> neighbours_target_type;

      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      
      std::auto_ptr<impl_type> tmp(pimpl);
      
      if (dynamic_cast<const neighbours_source_type*>(x.pimpl))
	pimpl = new neighbours_source_type();
      else
	pimpl = new neighbours_target_type();
      
      return *this;
    }
    
    template <typename FeaturePrefix, typename Feature>
    inline
    bool equal_prefix(const FeaturePrefix& prefix, const Feature& x)
    {
      return x.size() >= prefix.size() && std::equal(prefix.begin(), prefix.end(), x.begin());
    }

    
    void Neighbours::operator()(state_ptr_type& state,
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
      
      pimpl->neighbours_score(state, states, edge, features);
    }
    
    void Neighbours::operator()(const state_ptr_type& state,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates) const
    {
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();

      pimpl->neighbours_final_score(state, features);
    }

    void Neighbours::initialize()
    {
      pimpl->clear();
    }
  };
};
