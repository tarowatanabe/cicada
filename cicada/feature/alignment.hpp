// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__ALIGNMENT__HPP__
#define __CICADA__FEATURE__ALIGNMENT__HPP__ 1

#include <cicada/cluster.hpp>
#include <cicada/feature_function.hpp>

namespace cicada
{
  namespace feature
  {
    class RelativePosition : public FeatureFunction
    {
    public:
      typedef Cluster cluster_type;
      
    private:
      typedef FeatureFunction base_type;

    public:
      RelativePosition(const std::string& parameter);
      RelativePosition(const RelativePosition& x)
	: base_type(static_cast<const base_type&>(x)), cluster(0)
      {
	if (x.cluster)
	  cluster = &cluster_type::create(x.cluster->algorithm());
      }
      RelativePosition& operator=(const RelativePosition& x)
      {
	static_cast<base_type&>(*this) = static_cast<const base_type&>(x);
	
	cluster = 0;
	if (x.cluster)
	  cluster = &cluster_type::create(x.cluster->algorithm());
	return *this;
      }
      
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
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new RelativePosition(*this)); }
      
      virtual void assign(const size_type& id,
			  const hypergraph_type& hypergraph,
			  const lattice_type& lattice,
			  const span_set_type& spans,
			  const sentence_set_type& targets,
			  const ngram_count_set_type& ngram_counts);
      
    private:
      const cluster_type* cluster;
      const sentence_type* sentence;
      int source_length;
      int target_length;
    };
  };
};


#endif
