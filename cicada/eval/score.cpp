
#include "eval/score.hpp"
#include "eval/bleu.hpp"

#include "parameter.hpp"

namespace cicada
{
  namespace eval
  {

    std::string Scorer::lists()
    {
      return "bleu,order=<order, default=4>\n";
    }
    
    Scorer::scorer_ptr_type Scorer::create(const std::string& parameter)
    {
      typedef cicada::Parameter parameter_type;

      const parameter_type param(parameter);
      
      if (param.name() == "bleu" || param.name() == "bleu-linear") {
	int order = 4;
	parameter_type::const_iterator iter = param.find("order");
	if (iter != param.end())
	  order = boost::lexical_cast<int>(iter->second);
	
	return scorer_ptr_type(new BleuScorer(order));
      } else
	throw std::runtime_error("unknown scorer" + param.name());
      
    }
    
  };
};
