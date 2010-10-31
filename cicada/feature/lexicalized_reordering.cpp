
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

      LexicalizedReorderingImpl(const std::string& parameter) : model(parameter), lattice(0) {}
      LexicalizedReorderingImpl(const LexicalizedReorderingImpl& x) : model(x.model), lattice(0) {}
      LexicalizedReorderingImpl& operator=(const LexicalizedReorderingImpl& x)
      {
	model = x.model;
	return *this;
      }

      double reordering_score_next(const int prev_first, const int prev_last,
				   const int next_first, const int next_last,
				   const size_type node) const
      {
	const model_type::feature_set_type& features = model[node];

	if (! model.bidirectional)
	  return 0.0;
	else {
	  if (model.monotonicity) {
	    if (prev_last == next_first)
	      return features[2];
	    else
	      return features[3];
	  } else {
	    if (prev_last == next_first)
	      return features[3];
	    else if (next_last == prev_first)
	      return features[4];
	    else
	      return features[5];
	  }
	}
      }
      
      double reordering_score(const int prev_first, const int prev_last,
			      const int next_first, const int next_last,
			      const size_type node) const
      {
	const model_type::feature_set_type& features = model[node];
	
	if (model.monotonicity) {
	  if (prev_last == next_first)
	    return features[0];
	  else
	    return features[1];
	} else {
	  if (prev_last == next_first)
	    return features[0];
	  else if (next_last == prev_first)
	    return features[1];
	  else
	    return features[2];
	}
      }
      
      
      double lexicalized_reordering_score(state_ptr_type& state,
					 const state_ptr_set_type& states,
					 const edge_type& edge) const
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
	  
	  return reordering_score(0, 0, edge.first, edge.last, *node);
	} else if (states.size() == 1) {
	  // it is only for the goal state...
	  const int*       span_antecedent = reinterpret_cast<const int*>(states[0]);
	  const size_type* node_antecedent = reinterpret_cast<const size_type*>(span_antecedent + 2);

	  span[0] = span_antecedent[0];
	  span[1] = span_antecedent[1];
	  *node   = *node_antecedent;
	  return 0.0;
	} else if (states.size() == 2) {
	  const int*       span_antecedent = reinterpret_cast<const int*>(states[0]);
	  const size_type* node_antecedent = reinterpret_cast<const size_type*>(span_antecedent + 2);
	  
	  const int*       span_phrase     = reinterpret_cast<const int*>(states[1]);
	  const size_type* node_phrase     = reinterpret_cast<const size_type*>(span_phrase + 2);
	  
	  span[0] = span_phrase[0];
	  span[1] = span_phrase[1];
	  *node   = *node_phrase;
	  
	  // make adjustment, since this span-phrase is not a initial phrase..
	  const double score_adjust = reordering_score(0, 0,
						       span_phrase[0], span_phrase[1],
						       *node_phrase);
	  const double score        = reordering_score(span_antecedent[0], span_antecedent[1],
						       span_phrase[0], span_phrase[1],
						       *node_phrase);
	  const double score_next   = reordering_score_next(span_antecedent[0], span_antecedent[1],
							    span_phrase[0], span_phrase[1],
							    *node_antecedent);
	  
	  return score + score_next - score_adjust;
	} else
	  throw std::runtime_error("we do not support non-phrasal composed hypergraph");
      }

      double lexicalized_reordering_final_score(const state_ptr_type& state) const
      {
	const int* span = reinterpret_cast<const int*>(state);
	
	return (lattice ? - lattice->shortest_distance(span[1], lattice->size()) : 0);
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
      base_type::__feature_name = std::string("lexicalizedreordering");
      
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
      double score = pimpl->lexicalized_reordering_score(state, states, edge);
      
      if (final)
	score += pimpl->lexicalized_reordering_final_score(state);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
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
