//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/distortion.hpp"
#include "cicada/parameter.hpp"

#include "utils/indexed_set.hpp"
#include "utils/compact_trie_dense.hpp"
#include "utils/lexical_cast.hpp"

#include <boost/tuple/tuple.hpp>

namespace cicada
{
  namespace feature
  {
    
    class DistortionImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;
      typedef cicada::Lattice  lattice_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef attribute_set_type::attribute_type attribute_type;

      DistortionImpl()
	: lattice(0),
	  attr_phrase_span_first("phrase-span-first"),
	  attr_phrase_span_last("phrase-span-last") {}
      DistortionImpl(const DistortionImpl& x)
	: lattice(0),
	  attr_phrase_span_first("phrase-span-first"),
	  attr_phrase_span_last("phrase-span-last") {}
      DistortionImpl& operator=(const DistortionImpl& x)
      {
	return *this;
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
      
      double distortion_score(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge) const
      {
	int* span = reinterpret_cast<int*>(state);
	
	if (states.empty()) {
	  // How do we capture initial phrase....???
	  span[0] = phrase_span(edge.attributes, attr_phrase_span_first);
	  span[1] = phrase_span(edge.attributes, attr_phrase_span_last);
	  
	  return (lattice ? - lattice->shortest_distance(0, span[0]) : - span[0]);
	} else if (states.size() == 1) {
	  // it is only for the goal state...
	  const int* span_antecedent = reinterpret_cast<const int*>(states[0]);

	  span[0] = span_antecedent[0];
	  span[1] = span_antecedent[1];
	  return 0.0;
	} else if (states.size() == 2) {
	  const int* span_antecedent = reinterpret_cast<const int*>(states[0]);
	  const int* span_phrase     = reinterpret_cast<const int*>(states[1]);
	  
	  span[0] = span_phrase[0];
	  span[1] = span_phrase[1];
	  
	  // make adjustment, since this span-phrase is not a initial phrase..
	  const int score_adjust = (lattice ? lattice->shortest_distance(0, span[0]) : span[0]);

	  if (lattice) {
	    if (span_antecedent[1] == span_phrase[0])
	      return score_adjust;
	    else if (span_antecedent[1] < span_phrase[0])
	      return - lattice->shortest_distance(span_antecedent[1], span_phrase[0]) + score_adjust;
	    else
	      return - lattice->shortest_distance(span_phrase[0], span_antecedent[1]) + score_adjust;
	  } else
	    return - utils::bithack::abs(span_antecedent[1] - span_phrase[0]) + score_adjust;
	} else
	  throw std::runtime_error("we do not support non-phrasal composed hypergraph");
      }

      double distortion_final_score(const state_ptr_type& state) const
      {
	const int* span = reinterpret_cast<const int*>(state);
	
	return (lattice ? - lattice->shortest_distance(span[1], lattice->size()) : 0);
      }
      
      void assign(const lattice_type& __lattice)
      {
	lattice = &__lattice;
      }
      
      const lattice_type* lattice;
      
      const attribute_type attr_phrase_span_first;
      const attribute_type attr_phrase_span_last;
    };

    
    Distortion::Distortion(const std::string& parameter)
      : pimpl(0)
    {
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);

      if (param.name() != "distortion")
	throw std::runtime_error("is this really distortion feature function? " + parameter);
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter)
	std::cerr << "WARNING: unsupported parameter for distortion: " << piter->first << "=" << piter->second << std::endl;
      
      // distotion context: span = [first, last)
      base_type::__state_size = sizeof(int) * 2;
      base_type::__feature_name = std::string("distortion");
      
      pimpl = new impl_type();
    }
    
    Distortion::~Distortion() { std::auto_ptr<impl_type> tmp(pimpl); }

    
    Distortion::Distortion(const Distortion& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    { }
    
    
    Distortion& Distortion::operator=(const Distortion& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Distortion::apply(state_ptr_type& state,
			   const state_ptr_set_type& states,
			   const edge_type& edge,
			   feature_set_type& features,
			   feature_set_type& estimates,
			   const bool final) const
    {
      double score = pimpl->distortion_score(state, states, edge);
      
      if (final)
	score += pimpl->distortion_final_score(state);
      
      if (score != 0.0)
	features[base_type::feature_name()] = score;
      else
	features.erase(base_type::feature_name());
    }
    
    void Distortion::apply_coarse(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void Distortion::apply_predict(state_ptr_type& state,
				   const state_ptr_set_type& states,
				   const edge_type& edge,
				   feature_set_type& features,
				   feature_set_type& estimates,
				   const bool final) const
    {}
    
    void Distortion::apply_scan(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				const int dot,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {}
    void Distortion::apply_complete(state_ptr_type& state,
				    const state_ptr_set_type& states,
				    const edge_type& edge,
				    feature_set_type& features,
				    feature_set_type& estimates,
				    const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }

    
    void Distortion::assign(const size_type& id,
			    const hypergraph_type& hypergraph,
			    const lattice_type& lattice,
			    const span_set_type& spans,
			    const sentence_set_type& targets,
			    const ngram_count_set_type& ngram_counts)
    {
      pimpl->assign(lattice);
    }

    void Distortion::initialize()
    {
      
    }
  };
};
