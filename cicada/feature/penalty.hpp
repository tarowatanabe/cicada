// -*- mode: c++ -*-

#ifndef __CICADA__FEATURE__PENALTY__HPP__
#define __CICADA__FEATURE__PENALTY__HPP__ 1

namespace cicada
{

  namespace feature
  {
    class TargetWordPenalty : public FeatureFunction
    {
    public:
      TargetWordPenalty() : FeatureFunction(0, "target-word-penalty") { }
      
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
      
      void apply_estimate(const edge_type& edge,
			  feature_set_type& features) const
      {
	int count = 0;
	rule_type::symbol_set_type::const_iterator titer_end = edge.rule->target.end();
	for (rule_type::symbol_set_type::const_iterator titer = edge.rule->target.begin(); titer != titer_end; ++ titer)
	  count -= (titer->is_terminal() && *titer != vocab_type::EPSILON);
	
	if (count)
	  features[feature_name()] = count;
      }


      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new TargetWordPenalty(*this)); }
    };

    class SourceWordPenalty : public FeatureFunction
    {
    public:
      SourceWordPenalty() : FeatureFunction(0, "source-word-penalty") { }
      
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
      
      void apply_estimate(const edge_type& edge,
			  feature_set_type& features) const
      {
	int count = 0;
	rule_type::symbol_set_type::const_iterator siter_end = edge.rule->source.end();
	for (rule_type::symbol_set_type::const_iterator siter = edge.rule->source.begin(); siter != siter_end; ++ siter)
	  count -= (siter->is_terminal() && *siter != vocab_type::EPSILON);
	
	if (count)
	  features[feature_name()] = count;
      }

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new SourceWordPenalty(*this)); }
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
