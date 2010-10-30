
#include <utility>
#include <memory>

#include "cicada/feature/distortion.hpp"
#include "cicada/parameter.hpp"
#include "cicada/cluster.hpp"
#include "cicada/stemmer.hpp"

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

      typedef rule_type::symbol_set_type phrase_type;

      DistortionImpl() : lattice(0) {}
      DistortionImpl(const DistortionImpl& x) : lattice(0) {}
      DistortionImpl& operator=(const DistortionImpl& x)
      {
	return *this;
      }
      
      double distortion_score(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge) const
      {
	int* span = reinterpret_cast<int*>(state);
	
	if (states.empty()) {
	  span[0] = edge.first;
	  span[1] = edge.last;
	  return 0.0;
	} else {
	  
	}
      }

      double distortion_final_score(const state_ptr_type& state) const
      {
	
      }
      
      void assign(const lattice_type& __lattice)
      {
	lattice = &__lattice;
      }
      
      const lattice_type* lattice;
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
