//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/sgml_tag.hpp"
#include "cicada/parameter.hpp"
#include "cicada/weight_vector.hpp"

#include "utils/lexical_cast.hpp"
#include "utils/hashmurmur.hpp"
#include "utils/compress_stream.hpp"
#include "utils/piece.hpp"

namespace cicada
{
  namespace feature
  {
    
    class SGMLTagImpl
    {
    public:
      typedef cicada::Symbol   symbol_type;
      typedef cicada::Vocab    vocab_type;
      typedef cicada::Sentence sentence_type;

      typedef cicada::HyperGraph hypergraph_type;
      typedef cicada::Lattice    lattice_type;
      
      typedef cicada::FeatureFunction feature_function_type;
      
      typedef feature_function_type::state_ptr_type     state_ptr_type;
      typedef feature_function_type::state_ptr_set_type state_ptr_set_type;
      
      typedef feature_function_type::edge_type edge_type;

      typedef feature_function_type::feature_set_type   feature_set_type;
      typedef feature_function_type::attribute_set_type attribute_set_type;
      
      typedef attribute_set_type::attribute_type attribute_type;
      typedef feature_set_type::feature_type     feature_type;
      
      typedef feature_function_type::rule_type rule_type;

      typedef cicada::WeightVector<double, std::allocator<double> > weight_set_type;
      
      SGMLTagImpl() {}
      
      void sgml_tag_score(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features) const
      {
	
      }

      void estimate(const hypergraph_type& forest) const
      {
	
	
	
      }
      
    };
    
    SGMLTag::SGMLTag(const std::string& parameter)
      : pimpl(new impl_type())
    {
      typedef boost::filesystem::path path_type;
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "sgml-tag")
	throw std::runtime_error("is this really sgmltag feature function? " + parameter);

      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	std::cerr << "WARNING: unsupported parameter for sgml_tag: " << piter->first << "=" << piter->second << std::endl;
      }
      
      std::auto_ptr<impl_type> sgml_tag(new impl_type());
      
      
      pimpl = sgml_tag.release();
      
      base_type::__state_size = sizeof(impl_type::id_type) + sizeof(int) * 2;
      base_type::__feature_name = "sgml-tag";
    }
    
    SGMLTag::~SGMLTag() { std::auto_ptr<impl_type> tmp(pimpl); }

    SGMLTag::SGMLTag(const SGMLTag& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    SGMLTag& SGMLTag::operator=(const SGMLTag& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void SGMLTag::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      pimpl->sgml_tag_score(state, states, edge, features);
    }
    
    void SGMLTag::apply_coarse(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void SGMLTag::apply_predict(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void SGMLTag::apply_scan(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     const int dot,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void SGMLTag::apply_complete(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {}
    
    void SGMLTag::assign(const size_type& id,
			 const hypergraph_type& hypergraph,
			 const lattice_type& lattice,
			 const span_set_type& spans,
			 const sentence_set_type& targets,
			 const ngram_count_set_type& ngram_counts)
    {
      
    }
    
  };
};
