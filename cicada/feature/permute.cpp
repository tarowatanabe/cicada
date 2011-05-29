//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include <utility>
#include <memory>

#include "cicada/feature/permute.hpp"
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
    
    class PermuteImpl
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
      
      PermuteImpl() : feature("permute"), collapse(false) {}

      struct __rule_permute : public boost::static_visitor<bool>
      {
	bool operator()(const attribute_set_type::int_type& x) const { return x; }
	template <typename Tp>
	bool operator()(const Tp& x) const { return false; }
      };
      
      
      void permute_score(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features) const
      {
	if (collapse) {
	  attribute_set_type::const_iterator aiter_end = edge.attributes.end();
	  for (attribute_set_type::const_iterator aiter = edge.attributes.begin(); aiter != aiter_end; ++ aiter)
	    if (aiter->first.size() >= feature.size() && std::equal(feature.begin(), feature.end(), aiter->first.begin()))
	      if (boost::apply_visitor(__rule_permute(), aiter->second))
		features[feature] += weights[feature];
	} else {
	  attribute_set_type::const_iterator aiter_end = edge.attributes.end();
	  for (attribute_set_type::const_iterator aiter = edge.attributes.begin(); aiter != aiter_end; ++ aiter)
	    if (aiter->first.size() >= feature.size() && std::equal(feature.begin(), feature.end(), aiter->first.begin()))
	      if (boost::apply_visitor(__rule_permute(), aiter->second))
		features[static_cast<const std::string&>(aiter->first)] += 1.0;
	}
      }
      
      weight_set_type weights;
      feature_type    feature;
      bool            collapse;
    };
    
    Permute::Permute(const std::string& parameter)
      : pimpl(new impl_type())
    {
      typedef boost::filesystem::path path_type;
      typedef cicada::Parameter parameter_type;
      
      const parameter_type param(parameter);
      
      if (utils::ipiece(param.name()) != "permute")
	throw std::runtime_error("is this really permute feature function? " + parameter);

      bool      collapse = false;
      path_type path;
      
      for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	if (utils::ipiece(piter->first) == "collapse")
	  collapse = utils::lexical_cast<bool>(piter->second);
	else if (utils::ipiece(piter->first) == "weights")
	  path = piter->second;
	else
	  std::cerr << "WARNING: unsupported parameter for permute: " << piter->first << "=" << piter->second << std::endl;
      }

      std::auto_ptr<impl_type> permute(new impl_type());
      
      if (! path.empty()) {
	if (! boost::filesystem::exists(path))
	  throw std::runtime_error("no weight file? " + path.string());
	
	utils::compress_istream is(path, 1024 * 1024);
	is >> permute->weights;
      }
      
      permute->collapse = collapse;

      pimpl = permute.release();
      
      // if collapsed, this is not a sparse feature.. but since we always apply feature, this is non-sparse
      
      base_type::__state_size = 0;
      base_type::__feature_name = "permute";
      //base_type::__sparse_feature = ! collapse;
    }
    
    Permute::~Permute() { std::auto_ptr<impl_type> tmp(pimpl); }

    Permute::Permute(const Permute& x)
      : base_type(static_cast<const base_type&>(x)),
	pimpl(new impl_type(*x.pimpl))
    {}
    
    Permute& Permute::operator=(const Permute& x)
    {
      static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
      *pimpl = *x.pimpl;
      
      return *this;
    }
    
    void Permute::apply(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
    {
      features.erase_prefix(static_cast<const std::string&>(base_type::feature_name()));
      
      pimpl->permute_score(state, states, edge, features);
    }
    
    void Permute::apply_coarse(state_ptr_type& state,
			       const state_ptr_set_type& states,
			       const edge_type& edge,
			       feature_set_type& features,
			       feature_set_type& estimates,
			       const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void Permute::apply_predict(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const
    {
      apply(state, states, edge, features, estimates, final);
    }
    
    void Permute::apply_scan(state_ptr_type& state,
			     const state_ptr_set_type& states,
			     const edge_type& edge,
			     const int dot,
			     feature_set_type& features,
			     feature_set_type& estimates,
			     const bool final) const
    {}
    void Permute::apply_complete(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const
    {}    
  };
};
