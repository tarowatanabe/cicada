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
      
      void operator()(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      feature_set_type& features,
		      feature_set_type& estimates) const
      {
	int count = 0;
	rule_type::symbol_set_type::const_iterator titer_end = edge.rule->target.end();
	for (rule_type::symbol_set_type::const_iterator titer = edge.rule->target.begin(); titer != titer_end; ++ titer)
	  count -= (titer->is_terminal() && *titer != vocab_type::EPSILON);
	
	if (count)
	  features[feature_name()] = count;
      }

      void operator()(const state_ptr_type& state,
		      feature_set_type& features,
		      feature_set_type& estimates) const {}

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new TargetWordPenalty(*this)); }
    };

    class SourceWordPenalty : public FeatureFunction
    {
    public:
      SourceWordPenalty() : FeatureFunction(0, "source-word-penalty") { }
      
      void operator()(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      feature_set_type& features,
		      feature_set_type& estimates) const
      {
	int count = 0;
	rule_type::symbol_set_type::const_iterator titer_end = edge.rule->source.end();
	for (rule_type::symbol_set_type::const_iterator titer = edge.rule->source.begin(); titer != titer_end; ++ titer)
	  count -= (titer->is_terminal() && *titer != vocab_type::EPSILON);
	
	if (count)
	  features[feature_name()] = count;
      }

      void operator()(const state_ptr_type& state,
		      feature_set_type& features,
		      feature_set_type& estimates) const {}

      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new SourceWordPenalty(*this)); }
    };
    
    class RulePenalty : public FeatureFunction
    {
    public:
      RulePenalty() : FeatureFunction(0, "rule-penalty") { }
      
      void operator()(state_ptr_type& state,
		      const state_ptr_set_type& states,
		      const edge_type& edge,
		      feature_set_type& features,
		      feature_set_type& estimates) const
      {
	features[feature_name()] = -1;
      }
      
      void operator()(const state_ptr_type& state,
		      feature_set_type& features,
		      feature_set_type& estimates) const {}
      
      virtual feature_function_ptr_type clone() const { return feature_function_ptr_type(new RulePenalty(*this)); }
    };
    
  };
};

#endif
