// -*- mode: c++ -*-
//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FEATURE__PENALTY__HPP__
#define __CICADA__FEATURE__PENALTY__HPP__ 1

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

    private:      
      void apply_estimate(const edge_type& edge,
			  feature_set_type& features) const
      {
	int count = 0;
	rule_type::symbol_set_type::const_iterator titer_end = edge.rule->rhs.end();
	for (rule_type::symbol_set_type::const_iterator titer = edge.rule->rhs.begin(); titer != titer_end; ++ titer)
	  count -= (titer->is_terminal() && *titer != vocab_type::EPSILON);
	
	if (count)
	  features[feature_name()] = count;
      }


      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new WordPenalty(*this)); }
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
      
    private:
      void apply_estimate(const edge_type& edge,
			  feature_set_type& features) const
      {
	features[feature_name()] = -1;
      }
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new RulePenalty(*this)); }
    };
    
  };
};

#endif
