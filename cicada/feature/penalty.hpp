// -*- mode: c++ -*-
//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__PENALTY__HPP__
#define __CICADA__FEATURE__PENALTY__HPP__ 1

#include <vector>

#include <cicada/feature_function.hpp>

namespace cicada
{

  namespace feature
  {
    class WordPenalty : public FeatureFunction
    {
    public:
      WordPenalty() : FeatureFunction(0, "word-penalty") { }
      
      void apply(state_ptr_type& state,
		 const state_ptr_set_type& states,
		 const edge_type& edge,
		 feature_set_type& features,
		 feature_set_type& estimates,
		 const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_coarse(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
      {
	apply_estimate(edge, features);
      }

      void apply_predict(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
      {
	apply_estimate(edge, features);
      }

      void apply_scan(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      const int dot,
		      feature_set_type& features,
		      feature_set_type& estimates,
		      const bool final) const {}

      void apply_complete(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const {}
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new WordPenalty(*this)); }
      
    private:      
      void apply_estimate(const edge_type& edge,
			  feature_set_type& features) const
      {
	int count = 0;
	rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer)
	  count += (titer->is_terminal() && *titer != vocab_type::EPSILON && *titer != vocab_type::BOS && *titer != vocab_type::EOS);
	
	if (count)
	  features[feature_name()] = -count;
	else
	  features.erase(feature_name());
      }
    };

    
    class RulePenalty : public FeatureFunction
    {
    public:
      RulePenalty() : FeatureFunction(0, "rule-penalty") { }
      
      void apply(state_ptr_type& state,
		 const state_ptr_set_type& states,
		 const edge_type& edge,
		 feature_set_type& features,
		 feature_set_type& estimates,
		 const bool final) const
      {
	apply_estimate(edge, features);
      }
      void apply_coarse(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_predict(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_scan(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      const int dot,
		      feature_set_type& features,
		      feature_set_type& estimates,
		      const bool final) const {}
      
      void apply_complete(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const {}

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new RulePenalty(*this)); }
      
    private:
      void apply_estimate(const edge_type& edge,
			  feature_set_type& features) const
      {
	features[feature_name()] = -1;
      }
    };

    class ArityPenalty : public FeatureFunction
    {
    public:
      ArityPenalty() : FeatureFunction(0, "arity-penalty") { }
      
      void apply(state_ptr_type& state,
		 const state_ptr_set_type& states,
		 const edge_type& edge,
		 feature_set_type& features,
		 feature_set_type& estimates,
		 const bool final) const
      {
	apply_estimate(edge, features);
      }
      void apply_coarse(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_predict(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_scan(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      const int dot,
		      feature_set_type& features,
		      feature_set_type& estimates,
		      const bool final) const {}
      
      void apply_complete(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const {}

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new ArityPenalty(*this)); }
      
    private:
      void apply_estimate(const edge_type& edge, feature_set_type& features) const;
      
      typedef std::vector<feature_type, std::allocator<feature_type> > feature_name_set_type;
      
      feature_name_set_type names;
    };

    class GlueTreePenalty : public FeatureFunction
    {
    public:
      GlueTreePenalty() : FeatureFunction(0, "glue-tree-penalty"), attr_glue_tree("glue-tree") { }
      
      void apply(state_ptr_type& state,
		 const state_ptr_set_type& states,
		 const edge_type& edge,
		 feature_set_type& features,
		 feature_set_type& estimates,
		 const bool final) const
      {
	apply_estimate(edge, features);
      }
      void apply_coarse(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_predict(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_scan(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      const int dot,
		      feature_set_type& features,
		      feature_set_type& estimates,
		      const bool final) const {}
      
      void apply_complete(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const {}

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new GlueTreePenalty(*this)); }
      
    private:
      
      struct __glue_tree : public boost::static_visitor<bool>
      {
	bool operator()(const attribute_set_type::int_type& x) const { return x; }
	template <typename Tp>
	bool operator()(const Tp& x) const { return false; }
      };
      
      void apply_estimate(const edge_type& edge,
			  feature_set_type& features) const
      {
	attribute_set_type::const_iterator aiter = edge.attributes.find(attr_glue_tree);
	if (aiter != edge.attributes.end() && boost::apply_visitor(__glue_tree(), aiter->second))
	  features[feature_name()] = -1;
	else
	  features.erase(feature_name());
      }

      attribute_set_type::attribute_type attr_glue_tree;
    };

    class NonLatinPenalty : public FeatureFunction
    {
    public:
      NonLatinPenalty() : FeatureFunction(0, "non-latin-penalty") { }
      
      void apply(state_ptr_type& state,
		 const state_ptr_set_type& states,
		 const edge_type& edge,
		 feature_set_type& features,
		 feature_set_type& estimates,
		 const bool final) const
      {
	apply_estimate(edge, features);
      }
      void apply_coarse(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_predict(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_scan(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      const int dot,
		      feature_set_type& features,
		      feature_set_type& estimates,
		      const bool final) const {}
      
      void apply_complete(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const {}

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new NonLatinPenalty(*this)); }
      
    private:
      void apply_estimate(const edge_type& edge, feature_set_type& features) const;

    private:
      typedef std::vector<bool, std::allocator<bool> > non_latin_type;

      non_latin_type non_latin;
    };
    
    class InternalNodePenalty : public FeatureFunction
    {
    public:
      InternalNodePenalty() : FeatureFunction(0, "internal-node-penalty"),
			      attr_internal_node("internal-node") { }
      
      void apply(state_ptr_type& state,
		 const state_ptr_set_type& states,
		 const edge_type& edge,
		 feature_set_type& features,
		 feature_set_type& estimates,
		 const bool final) const
      {
	apply_estimate(edge, features);
      }
      void apply_coarse(state_ptr_type& state,
			const state_ptr_set_type& states,
			const edge_type& edge,
			feature_set_type& features,
			feature_set_type& estimates,
			const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_predict(state_ptr_type& state,
			 const state_ptr_set_type& states,
			 const edge_type& edge,
			 feature_set_type& features,
			 feature_set_type& estimates,
			 const bool final) const
      {
	apply_estimate(edge, features);
      }
      
      void apply_scan(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      const int dot,
		      feature_set_type& features,
		      feature_set_type& estimates,
		      const bool final) const {}
      
      void apply_complete(state_ptr_type& state,
			  const state_ptr_set_type& states,
			  const edge_type& edge,
			  feature_set_type& features,
			  feature_set_type& estimates,
			  const bool final) const {}

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new InternalNodePenalty(*this)); }
      
    private:
      void apply_estimate(const edge_type& edge, feature_set_type& features) const;
      
    private:
      typedef attribute_set_type::attribute_type attribute_type;
      
      attribute_type attr_internal_node;
    };

  };
};

#endif
