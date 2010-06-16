
#include "feature_function.hpp"
#include "parameter.hpp"

#include "feature/ngram.hpp"
#include "feature/penalty.hpp"

namespace cicada
{
  
  std::string FeatureFunction::lists()
  {
    static const char* desc = "\
ngram:file=<file>,order=<order>,name=feature-name(default: ngram)\n\
target-word-penalty\n\
soruce-word-penalty\n\
rule-penalty\n\
";
      return desc;
    
  }

  FeatureFunction::feature_function_ptr_type FeatureFunction::create(const std::string& parameter)
  {
    typedef cicada::Parameter parameter_type;
    
    const parameter_type param(parameter);
    
    if (param.name() == "ngram")
      return feature_function_ptr_type(new feature::NGram(parameter));
    else if (param.name() == "target-word-penalty")
      return feature_function_ptr_type(new feature::TargetWordPenalty());
    else if (param.name() == "source-word-penalty")
      return feature_function_ptr_type(new feature::SourceWordPenalty());
    else if (param.name() == "rule-penalty")
      return feature_function_ptr_type(new feature::RulePenalty());
    else
      throw std::runtime_error("unknown featuer: " + parameter);
    
  }
  
};
