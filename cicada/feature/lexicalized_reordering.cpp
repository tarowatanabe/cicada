
#include <utility>
#include <memory>

#include "cicada/feature/lexicalized_reordering.hpp"
#include "cicada/parameter.hpp"
#include "cicada/lexicalized_reordering.hpp"

#include "utils/indexed_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/lexical_cast.hpp"

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
      
      typedef cicada::LexicalizedReordering model_type;

      typedef model_type::size_type       size_type;
      typedef model_type::difference_type difference_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;
      typedef feature_function_type::feature_set_type feature_set_type;
      
      typedef feature_set_type::feature_type feature_type;
      typedef std::vector<feature_type, std::allocator<feature_type> > feature_list_type;
      
      LexicalizedReorderingImpl(const std::string& parameter)
	: model(parameter), feature_names(), lattice(0)
      {
	if (model.bidirectional) {
	  if (model.monotonicity) {
	    feature_names.push_back("lexicalized-reordering:forward-monotone");
	    feature_names.push_back("lexicalized-reordering:forward-others");
	    feature_names.push_back("lexicalized-reordering:backward-monotone");
	    feature_names.push_back("lexicalized-reordering:backward-others");
	  } else {
	    feature_names.push_back("lexicalized-reordering:forward-monotone");
	    feature_names.push_back("lexicalized-reordering:forward-swap");
	    feature_names.push_back("lexicalized-reordering:forward-discontinuous");
	    feature_names.push_back("lexicalized-reordering:backward-monotone");
	    feature_names.push_back("lexicalized-reordering:backward-swap");
	    feature_names.push_back("lexicalized-reordering:backward-discontinuous");
	  }
	} else {
	  if (model.monotonicity) {
	    feature_names.push_back("lexicalized-reordering:forward-monotone");
	    feature_names.push_back("lexicalized-reordering:forward-others");
	  } else {
	    feature_names.push_back("lexicalized-reordering:forward-monotone");
	    feature_names.push_back("lexicalized-reordering:forward-swap");
	    feature_names.push_back("lexicalized-reordering:forward-discontinuous");
	  }
	}
      }
      
      LexicalizedReorderingImpl(const LexicalizedReorderingImpl& x)
	: model(x.model), feature_names(x.feature_names), lattice(0) {}
      LexicalizedReorderingImpl& operator=(const LexicalizedReorderingImpl& x)
      {
	model = x.model;
	feature_names = x.feature_names;
	return *this;
      }

      void assign_feature(const feature_type& feature, const double& score, feature_set_type& features) const
      {
	if (score != 0.0)
	  features[feature] += score;
      }

      void reordering_score_next(const int prev_first, const int prev_last,
				 const int next_first, const int next_last,
				 const size_type node,
				 feature_set_type& features) const
      {
	if (model.bidirectional) {
	  if (model.monotonicity) {
	    if (prev_last == next_first)
	      assign_feature(feature_names[2], model[node][2], features);
	    else
	      assign_feature(feature_names[3], model[node][3], features);
	  } else {
	    if (prev_last == next_first)
	      assign_feature(feature_names[3], model[node][3], features);
	    else if (next_last == prev_first)
	      assign_feature(feature_names[4], model[node][4], features);
	    else
	      assign_feature(feature_names[5], model[node][5], features);
	  }
	}
      }
      
      void reordering_score(const int prev_first, const int prev_last,
			    const int next_first, const int next_last,
			    const size_type node,
			    feature_set_type& features) const
      {
	if (model.monotonicity) {
	  if (prev_last == next_first)
	    assign_feature(feature_names[0], model[node][0], features);
	  else
	    assign_feature(feature_names[1], model[node][1], features);
	} else {
	  if (prev_last == next_first)
	    assign_feature(feature_names[0], model[node][0], features);
	  else if (next_last == prev_first)
	    assign_feature(feature_names[1], model[node][1], features);
	  else
	    assign_feature(feature_names[2], model[node][2], features);
	}
      }
    
      void reordering_score_adjust(const int prev_first, const int prev_last,
				   const int next_first, const int next_last,
				   const size_type node,
				   feature_set_type& features) const
      {
	if (model.monotonicity) {
	  if (prev_last == next_first)
	    assign_feature(feature_names[0], - model[node][0], features);
	  else
	    assign_feature(feature_names[1], - model[node][1], features);
	} else {
	  if (prev_last == next_first)
	    assign_feature(feature_names[0], - model[node][0], features);
	  else if (next_last == prev_first)
	    assign_feature(feature_names[1], - model[node][1], features);
	  else
	    assign_feature(feature_names[2], - model[node][2], features);
	}
      }
      
      
      void lexicalized_reordering_score(state_ptr_type& state,
					const state_ptr_set_type& states,
					const edge_type& edge,
					feature_set_type& features) const
      {
	if (! lattice)
	  throw std::runtime_error("no input lattice?");

	int* span       = reinterpret_cast<int*>(state);
	size_type* node = reinterpret_cast<size_type*>(span + 2);
	
	if (states.empty()) {
	  // How do we capture initial phrase....???
	  span[0] = edge.first;
	  span[1] = edge.last;

	  sentence_type& phrase = const_cast<sentence_type&>(phrase_impl);
	  phrase.clear();
	  for (int pos = edge.first; pos != edge.last; ++ pos)
	    phrase.push_back(lattice->operator[](pos).front().label);
	  
	  *node = (model.fe
		   ? model.find(phrase.begin(), phrase.end(), edge.rule->rhs.begin(), edge.rule->rhs.end())
		   : model.find(phrase.begin(), phrase.end()));
	  
	  reordering_score(0, 0, edge.first, edge.last, *node, features);
	} else if (states.size() == 1) {
	  // it is only for the goal state...
	  const int*       span_antecedent = reinterpret_cast<const int*>(states[0]);
	  const size_type* node_antecedent = reinterpret_cast<const size_type*>(span_antecedent + 2);

	  span[0] = span_antecedent[0];
	  span[1] = span_antecedent[1];
	  *node   = *node_antecedent;
	} else if (states.size() == 2) {
	  const int*       span_antecedent = reinterpret_cast<const int*>(states[0]);
	  const size_type* node_antecedent = reinterpret_cast<const size_type*>(span_antecedent + 2);
	  
	  const int*       span_phrase     = reinterpret_cast<const int*>(states[1]);
	  const size_type* node_phrase     = reinterpret_cast<const size_type*>(span_phrase + 2);
	  
	  span[0] = span_phrase[0];
	  span[1] = span_phrase[1];
	  *node   = *node_phrase;
	  
	  // make adjustment, since this span-phrase is not a initial phrase..
	  reordering_score_adjust(0, 0, span_phrase[0], span_phrase[1], *node_phrase, features);
	  reordering_score(span_antecedent[0], span_antecedent[1], span_phrase[0], span_phrase[1], *node_phrase, features);
	  reordering_score_next(span_antecedent[0], span_antecedent[1], span_phrase[0], span_phrase[1], *node_antecedent, features);
	} else
	  throw std::runtime_error("we do not support non-phrasal composed hypergraph");
      }
      
      void lexicalized_reordering_final_score(const state_ptr_type& state,
						feature_set_type& features) const
      {
	const int* span       = reinterpret_cast<const int*>(state);
	const size_type* node = reinterpret_cast<const size_type*>(span + 2);
	
	reordering_score_next(span[0], span[1], lattice->size(), lattice->size(), *node, features);
      }
      
      void assign(const lattice_type& __lattice)
      {
	lattice = &__lattice;
	
	lattice_type::const_iterator liter_end = lattice->end();
	for (lattice_type::const_iterator liter = lattice->begin(); liter != liter_end; ++ liter)
	  if (liter->size() != 1)
	    throw std::runtime_error("we do not support non-linear lattice!");
      }
      
      model_type model;
      sentence_type phrase_impl;

      feature_list_type feature_names;
      
      const lattice_type* lattice;
    };

    
    LexicalizedReordering::LexicalizedReordering(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "lexicalized-reordering")
	throw std::runtime_error("is this really lexicalized reordering feature function? " + parameter);
      
      std::string file;

      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (strcasecmp(piter->first.c_str(), "file") == 0)
	  file = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for lexicalizedreordering: " << piter->first << "=" << piter->second << std::endl;
      }

      if (file.empty())
	throw std::runtime_error("no lexicalized reordering model?");
      
      // distotion context: span = [first, last)
      base_type::__state_size = sizeof(int) * 2  + sizeof(impl_type::size_type);
      base_type::__feature_name = std::string("lexicalized-reordering");
      
      pimpl = new impl_type(file);
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
      
    }
  };
};
