// -*- mode: c++ -*-
//
//  Copyright(C) 2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__DEPENDENCY__HPP__
#define __CICADA__FEATURE__DEPENDENCY__HPP__ 1

#include <cicada/cluster_stemmer.hpp>
#include <cicada/feature_function.hpp>

namespace cicada
{
  namespace feature
  {
    class DependencyImpl;

    class Dependency : public FeatureFunction
    {
    public:
      typedef feature_set_type::feature_type     feature_type;
      typedef attribute_set_type::attribute_type attribute_type;

    private:
      typedef FeatureFunction base_type;
      typedef DependencyImpl impl_type;
      
    public:
      Dependency(const std::string& parameter);
      Dependency(const Dependency& x);
      ~Dependency();
      
      Dependency& operator=(const Dependency& x);
      
    private:
      Dependency() {}

    public:
      virtual void apply(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const;
      virtual void apply_coarse(state_ptr_type& state,
				const state_ptr_set_type& states,
				const edge_type& edge,
				feature_set_type& features,
				feature_set_type& estimates,
				const bool final) const;
      virtual void apply_predict(state_ptr_type& state,
				 const state_ptr_set_type& states,
				 const edge_type& edge,
				 feature_set_type& features,
				 feature_set_type& estimates,
				 const bool final) const;
      virtual void apply_scan(state_ptr_type& state,
			      const state_ptr_set_type& states,
			      const edge_type& edge,
			      const int dot,
			      feature_set_type& features,
			      feature_set_type& estimates,
			      const bool final) const;
      virtual void apply_complete(state_ptr_type& state,
				  const state_ptr_set_type& states,
				  const edge_type& edge,
				  feature_set_type& features,
				  feature_set_type& estimates,
				  const bool final) const;

      virtual void assign(const size_type& id,
			  const hypergraph_type& hypergraph,
			  const lattice_type& lattice,
			  const span_set_type& spans,
			  const sentence_set_type& targets,
			  const ngram_count_set_type& ngram_counts);

      virtual void initialize();
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new Dependency(*this)); }
      
    private:
      impl_type* pimpl;
    };
  };
};

#endif
