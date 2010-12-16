//
//  Copyright(C) 2010 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "eval/score.hpp"
#include "eval/per.hpp"
#include "eval/wer.hpp"
#include "eval/sb.hpp"
#include "eval/sk.hpp"
#include "eval/ter.hpp"
#include "eval/bleu.hpp"
#include "eval/wlcs.hpp"
#include "eval/combined.hpp"

#include "stemmer.hpp"
#include "parameter.hpp"

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/lexical_cast.hpp"

#include <boost/filesystem.hpp>


namespace cicada
{
  namespace eval
  {

    const char* Scorer::lists()
    {
      static const char* desc = "\
combined: combined scorer\n\
\terror=[true|fase] error metric\n\
\treward=[true|fase] reward metric\n\
\tmetric=[scorer spec] i.e. metric=\"bleu:order=4\"\n\
\tweight=[weight for the scorer]\n\
bleu:\n\
\torder=<order, default=4> ngram order\n\
\ttokenizer=[tokenizer spec]\n\
per: position indenendent error rate\n\
\ttokenizer=[tokenizer spec]\n\
wer: word error rate\n\
\ttokenizer=[tokenizer spec]\n\
ter: translation error rate\n\
\ttokenizer=[tokenizer spec]\n\
\tmatch=approximated match cost (default 0.2)\n\
\tsubstitution=substitution cost (default 1)\n\
\tinsertion=insertion cost (default 1)\n\
\tdeletion=deletion cost (default 1)\n\
\tshift=shift cost (default 1)\n\
sk: string kernel\n\
\tp=order of string kernel (default 4)\n\
\tdecay=decay factor for string kernel (default 0.8)\n\
\ttokenizer=[tokenizer spec]\n\
wlcs: weighted longest common subsequence\n\
\talpha=length factor, k^alpha (default 1.0)\n\
\ttokenizer=[tokenizer spec]\n\
sb: skip bigram\n\
\twindow=window size (default 4, < 0 for infinity, == 0 for non-skip bigram)\n\
\ttokenizer=[tokenizer spec]\n\
";

      return desc;
    }
    
    Scorer::scorer_ptr_type Scorer::create(const std::string& parameter)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);

      scorer_ptr_type scorer;

      if (param.name() == "combined") {
	bool error = false;
	bool reward = false;

	std::vector<std::string, std::allocator<std::string> > metrics;
	std::vector<double, std::allocator<double> > weights;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "error") == 0)
	    error = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "reward") == 0)
	    reward = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "metric") == 0)
	    metrics.push_back(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "weight") == 0)
	    weights.push_back(boost::lexical_cast<double>(piter->second));
	  else
	    std::cerr << "WARNING: unsupported parameter for combined: " << piter->first << "=" << piter->second << std::endl;
	}

	if (int(error) + reward != 1)
	  throw std::runtime_error("specify whether this is error or reward metric");

	if (metrics.empty())
	  throw std::runtime_error("no metrics?");
	if (weights.empty())
	  throw std::runtime_error("no weights?");
	if (metrics.size() != weights.size())
	  throw std::runtime_error("metrics and weights size do not match");
	
	std::auto_ptr<CombinedScorer> combined(new CombinedScorer());
	combined->error = error;
	
	for (size_t i = 0; i != metrics.size(); ++ i) {
	  combined->scorers.push_back(create(metrics[i]));
	  combined->weights.push_back(weights[i]);
	}
	
	scorer.reset(combined.release());
	
      } else if (param.name() == "bleu" || param.name() == "bleu-linear") {
	int  order = 4;
	const tokenizer_type* tokenizer = 0;
	bool exact = false;

	bool yield_source = false;
	bool yield_target = false;
	
	std::string name;
	path_type   refset_file;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "order") == 0)
	    order = boost::lexical_cast<int>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "exact") == 0)
	    exact = utils::lexical_cast<bool>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "tokenizer") == 0)
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "name") == 0)
	    name = piter->second;
	  else if (strcasecmp(piter->first.c_str(), "refset") == 0)
	    refset_file = piter->second;
	  else if (strcasecmp(piter->first.c_str(), "yield") == 0) {
	    const std::string& yield = piter->second;
	    
	    if (strcasecmp(yield.c_str(), "source") == 0)
	      yield_source = true;
	    else if (strcasecmp(yield.c_str(), "target") == 0)
	      yield_target = true;
	    else
	      throw std::runtime_error("unknown parameter: " + parameter);
	    
	  } else
	    std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new BleuScorer(order));
	scorer->tokenizer = tokenizer;
      } else if (param.name() == "per") {
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "tokenizer") == 0)
	    tokenizer = &tokenizer_type::create(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for per: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new PERScorer());
	scorer->tokenizer = tokenizer;
      } else if (param.name() == "wer") {
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "tokenizer") == 0)
	    tokenizer = &tokenizer_type::create(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for wer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new WERScorer());
	scorer->tokenizer = tokenizer;
      } else if (param.name() == "ter") {
	const tokenizer_type* tokenizer = 0;

	double match = 0.2;
	double substitution = 1.0;
	double insertion = 1.0;
	double deletion = 1.0;
	double shift = 1.0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "tokenizer") == 0)
	    tokenizer = &tokenizer_type::create(piter->second);
	  if (strcasecmp(piter->first.c_str(), "match") == 0)
	    match = boost::lexical_cast<double>(piter->second);
	  if (strcasecmp(piter->first.c_str(), "substitution") == 0)
	    substitution = boost::lexical_cast<double>(piter->second);
	  if (strcasecmp(piter->first.c_str(), "insertion") == 0)
	    insertion = boost::lexical_cast<double>(piter->second);
	  if (strcasecmp(piter->first.c_str(), "deletion") == 0)
	    deletion = boost::lexical_cast<double>(piter->second);
	  if (strcasecmp(piter->first.c_str(), "shift") == 0)
	    shift = boost::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for ter: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new TERScorer(match, substitution, insertion, deletion, shift));
	scorer->tokenizer = tokenizer;
      } else if (param.name() == "sk") {
	int p = 4;
	double decay = 0.8;
	
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "tokenizer") == 0)
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "p") == 0)
	    p = boost::lexical_cast<int>(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "decay") == 0)
	    decay = boost::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for sk: " << piter->first << "=" << piter->second << std::endl;
	}

	scorer = scorer_ptr_type(new SKScorer(p, decay));
	scorer->tokenizer = tokenizer;
      } else if (param.name() == "wlcs") {
	double alpha = 1.0;
	
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "tokenizer") == 0)
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "alpha") == 0)
	    alpha = boost::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for wlcs: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new WLCSScorer(alpha));
	scorer->tokenizer = tokenizer;
      } else if (param.name() == "sb") {
	int window = 4;
	
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (strcasecmp(piter->first.c_str(), "tokenizer") == 0)
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (strcasecmp(piter->first.c_str(), "window") == 0)
	    window = boost::lexical_cast<int>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for sb: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new SBScorer(window));
	scorer->tokenizer = tokenizer;
      } else
	throw std::runtime_error("unknown scorer" + param.name());
      
      return scorer;
    }
    
  };
};
