//
//  Copyright(C) 2010-2011 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#include "eval/decode.hpp"

#include "eval/score.hpp"
#include "eval/per.hpp"
#include "eval/wer.hpp"
#include "eval/sb.hpp"
#include "eval/sk.hpp"
#include "eval/ter.hpp"
#include "eval/bleu.hpp"
#include "eval/wlcs.hpp"
#include "eval/combined.hpp"
#include "eval/parseval.hpp"

#include "stemmer.hpp"
#include "parameter.hpp"

#include "matcher.hpp"

#include <boost/thread.hpp>

#include "utils/config.hpp"
#include "utils/lexical_cast.hpp"
#include "utils/piece.hpp"

#include <boost/filesystem.hpp>


namespace cicada
{
  namespace eval
  {

    Score::score_ptr_type Score::decode(std::string::const_iterator& iter, std::string::const_iterator end)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      boost::spirit::qi::rule<std::string::const_iterator, std::string(), standard::space_type> quoted;
      
      quoted %= "\"" >> qi::lexeme[+(~standard::char_('\"'))] >> "\"";
      
      std::pair<std::string, std::string> scorer;

      std::string::const_iterator iter_saved = iter;
      
      const bool result = qi::phrase_parse(iter, end, (qi::lit('{') >> quoted >> qi::lit(':') >> quoted), standard::space, scorer);
      if (! result || scorer.first != "eval" || scorer.second.empty())
	return score_ptr_type();
      
      iter = iter_saved;
      
      if (scorer.second == "bleu")
	return Bleu::decode(iter, end);
      else if (scorer.second == "wer")
	return WER::decode(iter, end);
      else if (scorer.second == "per")
	return PER::decode(iter, end);
      else if (scorer.second == "ter")
	return TER::decode(iter, end);
      else if (scorer.second == "sk")
	return SK::decode(iter, end);
      else if (scorer.second == "sb")
	return SB::decode(iter, end);
      else if (scorer.second == "wlcs")
	return WLCS::decode(iter, end);
      else if (scorer.second == "combined")
	return Combined::decode(iter, end);
      else if (scorer.second == "parseval")
	return Parseval::decode(iter, end);
      else
	return score_ptr_type();
    }

    Score::score_ptr_type Score::decode(const utils::piece& encoded)
    {
      std::string::const_iterator iter(encoded.begin());
      std::string::const_iterator iter_end(encoded.end());
      
      return decode(iter, iter_end);
    }
    
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
\tskip-sgml-tag=[true|false] skip sgml tags\n\
per: position indenendent error rate\n\
\ttokenizer=[tokenizer spec]\n			\
wer: word error rate\n				\
\ttokenizer=[tokenizer spec]\n\
\tmatcher=[matcher spec] approximate matching\n\
\tmatch=approximated match cost (default 0.2)\n\
\tsubstitution=substitution cost (default 1)\n\
\tinsertion=insertion cost (default 1)\n\
\tdeletion=deletion cost (default 1)\n\
ter: translation error rate\n\
\ttokenizer=[tokenizer spec]\n\
\tmatcher=[matcher spec] approximate matching\n\
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
parseval: parse evaluation\n\
\tignore=[category] ignored category\n\
\ttokenizer=[tokenizer spec]\n\
";

      return desc;
    }
    
    Scorer::scorer_ptr_type Scorer::create(const utils::piece& parameter)
    {
      typedef cicada::Parameter parameter_type;
      typedef boost::filesystem::path path_type;
      
      const parameter_type param(parameter);

      scorer_ptr_type scorer;

      if (utils::ipiece(param.name()) == "combined") {
	bool error = false;
	bool reward = false;

	std::vector<std::string, std::allocator<std::string> > metrics;
	std::vector<double, std::allocator<double> > weights;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "error")
	    error = utils::lexical_cast<bool>(piter->second);
	  else if (utils::ipiece(piter->first) == "reward")
	    reward = utils::lexical_cast<bool>(piter->second);
	  else if (utils::ipiece(piter->first) == "metric")
	    metrics.push_back(piter->second);
	  else if (utils::ipiece(piter->first) == "weight")
	    weights.push_back(utils::lexical_cast<double>(piter->second));
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
	
      } else if (utils::ipiece(param.name()) == "bleu" || utils::ipiece(param.name()) == "bleu-linear") {
	int  order = 4;
	const tokenizer_type* tokenizer = 0;
	bool exact = false;
	bool skip_sgml_tag = false;

	bool yield_source = false;
	bool yield_target = false;
	
	std::string name;
	path_type   refset_file;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "order")
	    order = utils::lexical_cast<int>(piter->second);
	  else if (utils::ipiece(piter->first) == "exact")
	    exact = utils::lexical_cast<bool>(piter->second);
	  else if (utils::ipiece(piter->first) == "skip-sgml-tag")
	    skip_sgml_tag = utils::lexical_cast<bool>(piter->second);
	  else if (utils::ipiece(piter->first) == "tokenizer")
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (utils::ipiece(piter->first) == "name")
	    name = piter->second;
	  else if (utils::ipiece(piter->first) == "refset")
	    refset_file = piter->second;
	  else if (utils::ipiece(piter->first) == "yield") {
	    const utils::ipiece yield = piter->second;
	    
	    if (yield == "source")
	      yield_source = true;
	    else if (yield == "target")
	      yield_target = true;
	    else
	      throw std::runtime_error("unknown parameter: " + parameter);
	    
	  } else
	    std::cerr << "WARNING: unsupported parameter for bleu: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new BleuScorer(order));
	scorer->tokenizer = tokenizer;
      } else if (utils::ipiece(param.name()) == "per") {
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "tokenizer")
	    tokenizer = &tokenizer_type::create(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for per: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new PERScorer());
	scorer->tokenizer = tokenizer;
      } else if (utils::ipiece(param.name()) == "wer") {
	const tokenizer_type* tokenizer = 0;
	const Matcher* matcher = 0;
	
	WERScorer::weights_type weights;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "tokenizer")
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (utils::ipiece(piter->first) == "matcher")
	    matcher = &Matcher::create(piter->second);
	  else if (utils::ipiece(piter->first) == "match")
	    weights.match = utils::lexical_cast<double>(piter->second);
	  else if (utils::ipiece(piter->first) == "substitution")
	    weights.substitution = utils::lexical_cast<double>(piter->second);
	  else if (utils::ipiece(piter->first) == "insertion")
	    weights.insertion = utils::lexical_cast<double>(piter->second);
	  else if (utils::ipiece(piter->first) == "deletion")
	    weights.deletion = utils::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for wer: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new WERScorer(weights, matcher));
	scorer->tokenizer = tokenizer;
      } else if (utils::ipiece(param.name()) == "ter") {
	const tokenizer_type* tokenizer = 0;
	const Matcher* matcher = 0;
	
	TERScorer::weights_type weights;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "tokenizer")
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (utils::ipiece(piter->first) == "matcher")
	    matcher = &Matcher::create(piter->second);
	  else if (utils::ipiece(piter->first) == "match")
	    weights.match = utils::lexical_cast<double>(piter->second);
	  else if (utils::ipiece(piter->first) == "substitution")
	    weights.substitution = utils::lexical_cast<double>(piter->second);
	  else if (utils::ipiece(piter->first) == "insertion")
	    weights.insertion = utils::lexical_cast<double>(piter->second);
	  else if (utils::ipiece(piter->first) == "deletion")
	    weights.deletion = utils::lexical_cast<double>(piter->second);
	  else if (utils::ipiece(piter->first) == "shift")
	    weights.shift = utils::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for ter: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new TERScorer(weights, matcher));
	scorer->tokenizer = tokenizer;
      } else if (utils::ipiece(param.name()) == "sk") {
	int p = 4;
	double decay = 0.8;
	
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "tokenizer")
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (utils::ipiece(piter->first) == "p")
	    p = utils::lexical_cast<int>(piter->second);
	  else if (utils::ipiece(piter->first) == "decay")
	    decay = utils::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for sk: " << piter->first << "=" << piter->second << std::endl;
	}

	scorer = scorer_ptr_type(new SKScorer(p, decay));
	scorer->tokenizer = tokenizer;
      } else if (utils::ipiece(param.name()) == "wlcs") {
	double alpha = 1.0;
	
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "tokenizer")
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (utils::ipiece(piter->first) == "alpha")
	    alpha = utils::lexical_cast<double>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for wlcs: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new WLCSScorer(alpha));
	scorer->tokenizer = tokenizer;
      } else if (utils::ipiece(param.name()) == "sb") {
	int window = 4;
	
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "tokenizer")
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (utils::ipiece(piter->first) == "window")
	    window = utils::lexical_cast<int>(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for sb: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new SBScorer(window));
	scorer->tokenizer = tokenizer;
      } else if (utils::ipiece(param.name()) == "parseval") {
	std::vector<word_type, std::allocator<word_type> > ignored;
	
	const tokenizer_type* tokenizer = 0;
	
	for (parameter_type::const_iterator piter = param.begin(); piter != param.end(); ++ piter) {
	  if (utils::ipiece(piter->first) == "tokenizer")
	    tokenizer = &tokenizer_type::create(piter->second);
	  else if (utils::ipiece(piter->first) == "ignored")
	    ignored.push_back(piter->second);
	  else
	    std::cerr << "WARNING: unsupported parameter for parseval: " << piter->first << "=" << piter->second << std::endl;
	}
	
	scorer = scorer_ptr_type(new ParsevalScorer(ignored.begin(), ignored.end()));
	scorer->tokenizer = tokenizer;
      } else
	throw std::runtime_error("unknown scorer" + param.name());
      
      return scorer;
    }
    
  };
};
