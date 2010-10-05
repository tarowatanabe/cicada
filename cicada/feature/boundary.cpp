
#include <utility>
#include <memory>

#include "cicada/feature/boundary.hpp"
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
    
    class BoundaryImpl
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
      
      typedef std::pair<symbol_type, symbol_type> symbol_pair_type;
      typedef std::vector<symbol_pair_type, std::allocator<symbol_pair_type> > symbol_pair_set_type;

      typedef std::vector<int, std::allocator<int> > position_map_type;
      
      
      BoundaryImpl()
	: cluster_source(0), cluster_target(0),
	  stemmer_prefix_source(0), stemmer_prefix_target(0),
	  stemmer_suffix_source(0), stemmer_suffix_target(0),
	  stemmer_digits_source(0), stemmer_digits_target(0),
	  forced_feature(false) {}
      
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

      void compose_feature(std::string& name,
			   const std::string& source_prev, const std::string& source_next,
			   const std::string& target_prev, const std::string& target_next) const
      {
	name = "boundary:" + source_prev + '|' + source_next + '|' + target_prev + '|' + target_next;
      }
      
      void apply_feature(feature_set_type& features,
			 const symbol_type& source_prev, const symbol_type& source_next,
			 const symbol_type& target_prev, const symbol_type& target_next) const
      {
	if (cluster_source || cluster_target) {
	  const symbol_type source_prev_class = (cluster_source ? cluster_source->operator[](source_prev) : source_prev);
	  const symbol_type source_next_class = (cluster_source ? cluster_source->operator[](source_next) : source_next);
	  const symbol_type target_prev_class = (cluster_target ? cluster_target->operator[](target_prev) : target_prev);
	  const symbol_type target_next_class = (cluster_target ? cluster_target->operator[](target_next) : target_next);
	  
	  if (source_prev_class != source_prev || source_next_class != source_next
	      || target_prev_class != target_prev || target_next_class != target_next) {
	    
	    std::string name;
	    
	    compose_feature(name, source_prev_class, source_next_class, target_prev_class, target_next_class);
	    
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}
	
	if (stemmer_prefix_source || stemmer_prefix_target) {
	  const symbol_type source_prev_stemmed = (stemmer_prefix_source ? stemmer_prefix_source->operator[](source_prev) : source_prev);
	  const symbol_type source_next_stemmed = (stemmer_prefix_source ? stemmer_prefix_source->operator[](source_next) : source_next);
	  const symbol_type target_prev_stemmed = (stemmer_prefix_target ? stemmer_prefix_target->operator[](target_prev) : target_prev);
	  const symbol_type target_next_stemmed = (stemmer_prefix_target ? stemmer_prefix_target->operator[](target_next) : target_next);
	  
	  if (source_prev_stemmed != source_prev || source_next_stemmed != source_next
	      || target_prev_stemmed != target_prev || target_next_stemmed != target_next) {
	    
	    std::string name;
	    
	    compose_feature(name, source_prev_stemmed, source_next_stemmed, target_prev_stemmed, target_next_stemmed);
	    
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}
	
	if (stemmer_suffix_source || stemmer_suffix_target) {
	  const symbol_type source_prev_stemmed = (stemmer_suffix_source ? stemmer_suffix_source->operator[](source_prev) : source_prev);
	  const symbol_type source_next_stemmed = (stemmer_suffix_source ? stemmer_suffix_source->operator[](source_next) : source_next);
	  const symbol_type target_prev_stemmed = (stemmer_suffix_target ? stemmer_suffix_target->operator[](target_prev) : target_prev);
	  const symbol_type target_next_stemmed = (stemmer_suffix_target ? stemmer_suffix_target->operator[](target_next) : target_next);
	  
	  if (source_prev_stemmed != source_prev || source_next_stemmed != source_next
	      || target_prev_stemmed != target_prev || target_next_stemmed != target_next) {
	    
	    std::string name;
	    
	    compose_feature(name, source_prev_stemmed, source_next_stemmed, target_prev_stemmed, target_next_stemmed);
	    
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}
		
	if (stemmer_digits_source || stemmer_digits_target) {
	  const symbol_type source_prev_stemmed = (stemmer_digits_source ? stemmer_digits_source->operator[](source_prev) : source_prev);
	  const symbol_type source_next_stemmed = (stemmer_digits_source ? stemmer_digits_source->operator[](source_next) : source_next);
	  const symbol_type target_prev_stemmed = (stemmer_digits_target ? stemmer_digits_target->operator[](target_prev) : target_prev);
	  const symbol_type target_next_stemmed = (stemmer_digits_target ? stemmer_digits_target->operator[](target_next) : target_next);
	  
	  if (source_prev_stemmed != source_prev || source_next_stemmed != source_next
	      || target_prev_stemmed != target_prev || target_next_stemmed != target_next) {
	    
	    std::string name;
	    
	    compose_feature(name, source_prev_stemmed, source_next_stemmed, target_prev_stemmed, target_next_stemmed);
	    
	    if (forced_feature || feature_set_type::feature_type::exists(name))
	      features[name] += 1.0;
	  }
	}


	std::string name;
	
	compose_feature(name, source_prev, source_next, target_prev, target_next);
	
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

      cluster_type* cluster_source;
      cluster_type* cluster_target;

      stemmer_type* stemmer_prefix_source;
      stemmer_type* stemmer_prefix_target;
      stemmer_type* stemmer_suffix_source;
      stemmer_type* stemmer_suffix_target;
      stemmer_type* stemmer_digits_source;
      stemmer_type* stemmer_digits_target;
      
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

      int stemmer_prefix_source_size = 0;
      int stemmer_prefix_target_size = 0;
      int stemmer_suffix_source_size = 0;
      int stemmer_suffix_target_size = 0;
      
      bool stemmer_digits_source = false;
      bool stemmer_digits_target = false;
      
      boost::filesystem::path cluster_path_source;
      boost::filesystem::path cluster_path_target;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "cluster-source") == 0)
	  cluster_path_source = piter->second;
	else if (strcasecmp(piter->first.c_str(), "cluster-target") == 0)
	  cluster_path_target = piter->second;
	else if (strcasecmp(piter->first.c_str(), "prefix-source") == 0)
	  stemmer_prefix_source_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "prefix-target") == 0)
	  stemmer_prefix_target_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "suffix-source") == 0)
	  stemmer_suffix_source_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "suffix-target") == 0)
	  stemmer_suffix_target_size = boost::lexical_cast<int>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "digits-source") == 0)
	  stemmer_digits_source = utils::lexical_cast<bool>(piter->second);
	else if (strcasecmp(piter->first.c_str(), "digits-target") == 0)
	  stemmer_digits_target = utils::lexical_cast<bool>(piter->second);
	else
	  std::cerr << "WARNING: unsupported parameter for boundary: " << piter->first << "=" << piter->second << std::endl;
      }

      if (stemmer_prefix_source_size < 0)
	throw std::runtime_error("negative source prefix size?");
      if (stemmer_prefix_target_size < 0)
	throw std::runtime_error("negative target prefix size?");
      if (stemmer_suffix_source_size < 0)
	throw std::runtime_error("negative source suffix size?");
      if (stemmer_suffix_target_size < 0)
	throw std::runtime_error("negative target suffix size?");
      
      
      if (! cluster_path_source.empty()) {
	if (! boost::filesystem::exists(cluster_path_source))
	  throw std::runtime_error("no source cluster file: " + cluster_path_source.file_string());
	
	pimpl->cluster_source = &cicada::Cluster::create(cluster_path_source);
      }
      
      if (! cluster_path_target.empty()) {
	if (! boost::filesystem::exists(cluster_path_target))
	  throw std::runtime_error("no target cluster file: " + cluster_path_target.file_string());
	
	pimpl->cluster_target = &cicada::Cluster::create(cluster_path_target);
      }
      
      if (stemmer_prefix_source_size > 0)
	pimpl->stemmer_prefix_source = &cicada::Stemmer::create("prefix:size=" + boost::lexical_cast<std::string>(stemmer_prefix_source_size));
      if (stemmer_prefix_target_size > 0)
	pimpl->stemmer_prefix_target = &cicada::Stemmer::create("prefix:size=" + boost::lexical_cast<std::string>(stemmer_prefix_target_size));
      if (stemmer_suffix_source_size > 0)
	pimpl->stemmer_suffix_source = &cicada::Stemmer::create("suffix:size=" + boost::lexical_cast<std::string>(stemmer_suffix_source_size));
      if (stemmer_suffix_target_size > 0)
	pimpl->stemmer_suffix_target = &cicada::Stemmer::create("suffix:size=" + boost::lexical_cast<std::string>(stemmer_suffix_target_size));
      
      if (stemmer_digits_source)
	pimpl->stemmer_digits_source = &cicada::Stemmer::create("digits");
      if (stemmer_digits_target)
	pimpl->stemmer_digits_target = &cicada::Stemmer::create("digits");
      
      
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
    
    void Boundary::apply(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
    {
      const std::string& __feature_prefix = base_type::feature_name();
      for (feature_set_type::iterator fiter = features.begin(); fiter != features.end(); /**/)
	if (equal_prefix(__feature_prefix, fiter->first))
	  features.erase(fiter ++);
	else
	  ++ fiter;
      
      const_cast<impl_type*>(pimpl)->forced_feature = base_type::apply_feature();
      
      pimpl->boundary_score(state, states, edge, features);
      
      if (final)
	pimpl->boundary_final_score(state, features);
    }

    void Boundary::apply_coarse(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {

    }

    
    void Boundary::apply_predict(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {}
    void Boundary::apply_scan(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      const int dot,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const
    {}
    void Boundary::apply_complete(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }


    void Boundary::initialize()
    {
      
    }
  };
};
