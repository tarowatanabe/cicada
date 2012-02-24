//
//  Copyright(C) 2010-2012 Taro Watanabe <taro.watanabe@nict.go.jp>
//

#ifndef __CICADA__FILTER_EXTRACT_IMPL__HPP__
#define __CICADA__FILTER_EXTRACT_IMPL__HPP__ 1

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>
#include <boost/filesystem.hpp>
#include <boost/math/special_functions/expm1.hpp>

#include <string>
#include <vector>
#include <iostream>

#include <utils/hashmurmur.hpp>
#include <utils/alloc_vector.hpp>
#include <utils/compress_stream.hpp>
#include <utils/mathop.hpp>
#include <utils/dense_hash_map.hpp>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sentence.hpp>
#include <cicada/tree_rule.hpp>

struct RootCount
{
  typedef std::string label_type;
  typedef double count_type;
  typedef std::vector<count_type, std::allocator<count_type> > counts_type;
  
  label_type  label;
  counts_type counts;
  
  double observed_joint;
  double observed;
  
  RootCount() : label(), counts() {}
  RootCount(const label_type& __label) : label(__label), counts() {}
  
  void clear()
  {
    label.clear();
    counts.clear();
    observed_joint = 0;
    observed = 0;
  }
  
  friend
  size_t  hash_value(RootCount const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    return hasher_type()(x.label.begin(), x.label.end(), 0);
  }
  
  friend
  bool operator==(const RootCount& x, const RootCount& y) 
  {
    return x.label == y.label;
  }
  
  friend
  bool operator!=(const RootCount& x, const RootCount& y) 
  {
    return x.label != y.label;
  }
  
  friend
  bool operator<(const RootCount& x, const RootCount& y)
  {
    return x.label < y.label;
  }

  friend
  bool operator>(const RootCount& x, const RootCount& y)
  {
    return y < x;
  }
};

struct PhrasePair
{
  typedef std::string phrase_type;
  typedef double count_type;
  typedef std::vector<count_type, std::allocator<count_type> > counts_type;
  
  phrase_type source;
  phrase_type target;
  counts_type counts;
  counts_type counts_source;
  counts_type counts_target;
  double observed_source;
  double observed_target;
  double lexicon_target_source;
  double lexicon_source_target;
  
  PhrasePair()
    : source(), target(), counts(), counts_source(), counts_target() {}

  void clear()
  {
    source.clear();
    target.clear();
    counts.clear();
    counts_source.clear();
    counts_target.clear();
    observed_source = 0;
    observed_target = 0;
    lexicon_target_source = 0;
    lexicon_source_target = 0;
  }

  friend
  size_t  hash_value(PhrasePair const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;
   
    return hasher_type()(x.source.begin(), x.source.end(), hasher_type()(x.target.begin(), x.target.end(), 0));
  }


  friend
  bool operator==(const PhrasePair& x, const PhrasePair& y) 
  {
    return x.source == y.source && x.target == y.target;
  }
  
  friend
  bool operator!=(const PhrasePair& x, const PhrasePair& y) 
  {
    return x.source != y.source || x.target != y.target;
  }
  
  friend
  bool operator<(const PhrasePair& x, const PhrasePair& y)
  {
    return (x.source < y.source || (!(y.source < x.source) && x.target < y.target));
  }

  friend
  bool operator>(const PhrasePair& x, const PhrasePair& y)
  {
    return y < x;
  }

};

BOOST_FUSION_ADAPT_STRUCT(
			  RootCount,
			  (RootCount::label_type,  label)
			  (RootCount::counts_type, counts)
			  (RootCount::count_type,  observed_joint)
			  (RootCount::count_type, observed)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  PhrasePair,
			  (PhrasePair::phrase_type, source)
			  (PhrasePair::phrase_type, target)
			  (PhrasePair::counts_type, counts)
			  (PhrasePair::counts_type, counts_source)
			  (PhrasePair::counts_type, counts_target)
			  (PhrasePair::count_type,  observed_source)
			  (PhrasePair::count_type,  observed_target)
			  (PhrasePair::count_type,  lexicon_target_source)
			  (PhrasePair::count_type,  lexicon_source_target)
			  )

struct RootCountParser
{
  typedef RootCount root_count_type;
  
  typedef root_count_type::label_type  label_type;
  typedef root_count_type::counts_type counts_type;

  RootCountParser() : grammar() {}
  RootCountParser(const RootCountParser& x) : grammar() {}
  
  template <typename Iterator>
  struct root_count_parser : boost::spirit::qi::grammar<Iterator, root_count_type(), boost::spirit::standard::space_type>
  {
    root_count_parser() : root_count_parser::base_type(root_count)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      label %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];
      counts %= +qi::double_;
      root_count %= label >> "|||" >> counts >> "|||" >> qi::double_ >> qi::double_;
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> label;
    boost::spirit::qi::rule<Iterator, counts_type(), boost::spirit::standard::space_type> counts;
    boost::spirit::qi::rule<Iterator, root_count_type(), boost::spirit::standard::space_type> root_count;
  };

  bool operator()(std::istream& is, root_count_type& root_count)
  {
    root_count.clear();
    
    std::string line;
    if (! getline(is, line)) return false;
    
    return operator()(line, root_count);
  }
  
  bool operator()(const std::string& line, root_count_type& root_count)
  {
    root_count.clear();
    
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    const bool result =  boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, root_count);
    return result && iter == end;
  }
  
  root_count_parser<std::string::const_iterator> grammar;
};

struct PhrasePairParser
{
  typedef PhrasePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;
  
  PhrasePairParser() : grammar() {}
  PhrasePairParser(const PhrasePairParser& x) : grammar() {}

  template <typename Iterator>
  struct phrase_pair_parser : boost::spirit::qi::grammar<Iterator, phrase_pair_type(), boost::spirit::standard::space_type>
  {
    phrase_pair_parser() : phrase_pair_parser::base_type(phrase_pair)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      phrase %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];
      counts %= +qi::double_;
      
      phrase_pair %= (phrase
		      >> "|||" >> phrase
		      >> "|||" >> counts
		      >> "|||" >> counts
		      >> "|||" >> counts
		      >> "|||" >> qi::double_ >> qi::double_
		      >> "|||" >> qi::double_ >> qi::double_);
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, counts_type(), boost::spirit::standard::space_type> counts;
    boost::spirit::qi::rule<Iterator, phrase_pair_type(), boost::spirit::standard::space_type> phrase_pair;
  };

  bool operator()(std::istream& is, phrase_pair_type& phrase_pair)
  {
    phrase_pair.clear();
    
    std::string line;
    if (! getline(is, line)) return false;
    
    return operator()(line, phrase_pair);
  }
  
  bool operator()(const std::string& line, phrase_pair_type& phrase_pair)
  {
    phrase_pair.clear();
    
    std::string::const_iterator iter = line.begin();
    std::string::const_iterator end = line.end();
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, phrase_pair);
    
    if (result && iter == end)
      return true;
    else {
      std::cerr << "WARNING: parsing failed: " << line << std::endl;
      return false;
    }
  }
  
  phrase_pair_parser<std::string::const_iterator> grammar;
};

struct PhrasePairGenerator
{
  typedef PhrasePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;
  
  PhrasePairGenerator() : grammar() {}
  PhrasePairGenerator(const PhrasePairGenerator& x) : grammar() {}

  
  template <typename Iterator>
  struct phrase_pair_generator : boost::spirit::karma::grammar<Iterator, phrase_pair_type()>
  {
    phrase_pair_generator() : phrase_pair_generator::base_type(phrase_pair)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      
      counts %= double20 % ' ';
      phrase_pair %= (standard::string
		      << " ||| " << standard::string
		      << " ||| " << counts
		      << " ||| " << counts
		      << " ||| " << counts
		      << " ||| " << double20 << ' ' << double20
		      << " ||| " << double20 << ' ' << double20);
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 20;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double20;
    
    boost::spirit::karma::rule<Iterator, counts_type()> counts;
    boost::spirit::karma::rule<Iterator, phrase_pair_type()> phrase_pair;
  };

  typedef std::ostream_iterator<char> iterator_type;
  
  std::ostream& operator()(std::ostream& os, const phrase_pair_type& phrase_pair)
  {
    iterator_type iter(os);
    
    if (! boost::spirit::karma::generate(iter, grammar, phrase_pair))
      throw std::runtime_error("failed generation!");
    
    return os;
  }
  
  phrase_pair_generator<iterator_type> grammar;
};




struct ExtractRootPhrase
{
  const std::string& operator()(const std::string& phrase) const
  {
    static const std::string __label("[x]");
    return __label;
  }
};

struct ExtractRootSCFG
{
  // extract the first word...
  std::string operator()(const std::string& phrase) const
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();

    std::string label;
    
    const bool result = qi::phrase_parse(iter, end,
					 qi::lexeme[+(standard::char_ - standard::space)],
					 standard::space,
					 label);
    if (! result)
      throw std::runtime_error("no label?");
    
    return label;
  }
};

struct ExtractPhraseSCFG
{
  // extract the non-first word...
  std::pair<std::string, std::string> operator()(const std::string& phrase) const
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();
    
    std::pair<std::string, std::string> label_pair;
    
    const bool result = qi::phrase_parse(iter, end,
					 qi::lexeme[+(standard::char_ - standard::space)] >> qi::lexeme[+(standard::char_ - standard::space) >> *standard::char_],
					 standard::space,
					 label_pair);
    if (! result || iter != end)
      throw std::runtime_error("no label?");
    
    return label_pair;
  }
};

struct ExtractRootGHKM
{
  // extract the first word...
  std::string operator()(const std::string& phrase) const
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();
    
    std::string label;
    
    const bool result = qi::phrase_parse(iter, end,
					 qi::lexeme[+((standard::char_ - standard::space - '(') | "\\(")],
					 standard::space,
					 label);
    if (! result)
      throw std::runtime_error("no label?");
    
    return label;
  }
};

struct LexiconModel
{
  typedef cicada::Symbol word_type;
  typedef cicada::Vocab vocab_type;

  typedef boost::filesystem::path path_type;

  struct table_type
  {
    typedef google::dense_hash_map<word_type, double, boost::hash<word_type> , std::equal_to<word_type> > __table_type;

    table_type() : table() { table.set_empty_key(word_type()); }

    __table_type table;
  };

  typedef utils::alloc_vector<table_type, std::allocator<table_type> > table_set_type;

  LexiconModel() : smooth(1e-7), tables() {}
  LexiconModel(const path_type& path) : smooth(), tables() { open(path); }

  void open(const path_type& path)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    typedef boost::fusion::tuple<std::string, std::string, double> parsed_type;

    qi::rule<std::string::const_iterator, std::string(), standard::space_type> word;
    qi::rule<std::string::const_iterator, parsed_type(), standard::space_type> lexicon;

    word %= qi::lexeme[+(standard::char_ - standard::space)];
    lexicon %= word >> word >> qi::double_;

    smooth = std::numeric_limits<double>::infinity();
    tables.clear();

    utils::compress_istream is(path, 1024 * 1024);
    std::string line;
    parsed_type parsed;

    while (std::getline(is, line)) {
      std::string::const_iterator iter = line.begin();
      std::string::const_iterator end = line.end();

      boost::fusion::get<0>(parsed).clear();
      boost::fusion::get<1>(parsed).clear();

      const bool result = qi::phrase_parse(iter, end, lexicon, standard::space, parsed);
      if (! result || iter != end) continue;


      tables[word_type(boost::fusion::get<1>(parsed)).id()].table[boost::fusion::get<0>(parsed)] = boost::fusion::get<2>(parsed);

      smooth = std::min(smooth, boost::fusion::get<2>(parsed));
    }

  }

  double operator()(const word_type& source, const word_type& target) const
  {
    if (tables.exists(source.id())) {
      const table_type& table = tables[source.id()];

      table_type::__table_type::const_iterator iter = table.table.find(target);
      return (iter != table.table.end() ? iter->second : smooth);
    } else 
      return smooth;
  }

  double smooth;
  table_set_type tables;
};

struct LexiconBase
{
  typedef cicada::Sentence sentence_type;

  typedef LexiconModel lexicon_model_type;

  LexiconBase(const lexicon_model_type& __lexicon_source_target,
	      const lexicon_model_type& __lexicon_target_source)
    : lexicon_source_target(__lexicon_source_target),
      lexicon_target_source(__lexicon_target_source) {}

  const lexicon_model_type& lexicon_source_target;
  const lexicon_model_type& lexicon_target_source;
  
  std::pair<double, double> noisy_or() const
  {
    const size_t source_size = source.size();
    const size_t target_size = target.size();
    
    double score_source_target = 0.0;
    
    for (size_t trg = 0; trg != target_size; ++ trg)
      if (target[trg].is_terminal()) {
	double score = 0.0;
	for (size_t src = 0; src != source_size; ++ src)
	  if (source[src].is_terminal())
	    score += utils::mathop::log(1.0 - lexicon_source_target(source[src], target[trg]));
	
	//
	// 1.0 - exp(score) == - expm1(score)
	//
	score_source_target += utils::mathop::log(std::max(- boost::math::expm1(score), lexicon_source_target.smooth));
      }
    

    double score_target_source = 0.0;
    for (size_t src = 0; src != source_size; ++ src)
      if (source[src].is_terminal()) {
	double score = 0.0;
	for (size_t trg = 0; trg != target_size; ++ trg)
	  if (target[trg].is_terminal())
	    score += utils::mathop::log(1.0 - lexicon_target_source(target[trg], source[src]));
	
	score_target_source += utils::mathop::log(std::max(- boost::math::expm1(score), lexicon_target_source.smooth));
      }
    
    return std::make_pair(score_source_target, score_target_source);
  }
  
  std::pair<double, double> model1() const
  {
    const size_t source_size = source.size();
    const size_t target_size = target.size();
    
    double score_source_target = 0.0;
    
    for (size_t trg = 0; trg != target_size; ++ trg)
      if (target[trg].is_terminal()) {
	
	double score = lexicon_source_target(lexicon_model_type::vocab_type::EPSILON, target[trg]);
	for (size_t src = 0; src != source_size; ++ src)
	  if (source[src].is_terminal())
	    score += lexicon_source_target(source[src], target[trg]);
	
	score_source_target += utils::mathop::log(score);
      }
    
    double score_target_source = 0.0;
    
    for (size_t src = 0; src != source_size; ++ src)
      if (source[src].is_terminal()) {
	
	double score = lexicon_target_source(lexicon_model_type::vocab_type::EPSILON, source[src]);
	for (size_t trg = 0; trg != target_size; ++ trg)
	  if (target[trg].is_terminal())
	    score += lexicon_target_source(target[trg], source[src]);
	
	score_target_source += utils::mathop::log(score);
      }
    
    return std::make_pair(score_source_target, score_target_source);
  }
  
  std::pair<double, double> insertion_deletion(const double threshold_insertion=0.01,
					       const double threshold_deletion=0.01) const
  {
    const size_t source_size = source.size();
    const size_t target_size = target.size();
    
    double score_source_target = 0.0;
    
    for (size_t trg = 0; trg != target_size; ++ trg)
      if (target[trg].is_terminal()) {
	
	double score = 0.0;
	for (size_t src = 0; src != source_size; ++ src)
	  if (source[src].is_terminal())
	    score = std::max(score, lexicon_source_target(source[src], target[trg]));
	
	score_source_target -= (score < threshold_insertion);
      }
    
    double score_target_source = 0.0;
    
    for (size_t src = 0; src != source_size; ++ src)
      if (source[src].is_terminal()) {
	
	double score = 0.0;
	for (size_t trg = 0; trg != target_size; ++ trg)
	  if (target[trg].is_terminal())
	    score = std::max(score, lexicon_target_source(target[trg], source[src]));
	
	score_target_source -= (score < threshold_deletion);
      }
    
    return std::make_pair(score_source_target, score_target_source);
  }
  
  
  sentence_type source;
  sentence_type target;
};

struct LexiconPhrase : public LexiconBase
{

  LexiconPhrase(const lexicon_model_type& __lexicon_source_target,
		const lexicon_model_type& __lexicon_target_source)
    : LexiconBase(__lexicon_source_target, __lexicon_target_source) {}

  void assign_source(const std::string& phrase)
  {
    source.assign(phrase);
  }

  void assign_target(const std::string& phrase)
  {
    target.assign(phrase);
  }

};

struct LexiconSCFG : public LexiconBase
{

  LexiconSCFG(const lexicon_model_type& __lexicon_source_target,
	      const lexicon_model_type& __lexicon_target_source)
    : LexiconBase(__lexicon_source_target, __lexicon_target_source) {}

  void assign_source(const std::string& phrase)
  {
    source.assign(phrase);
    source.erase(source.begin());
  }

  void assign_target(const std::string& phrase)
  {
    target.assign(phrase);
    target.erase(target.begin());
  }
};

struct LexiconGHKM : public LexiconBase
{

  LexiconGHKM(const lexicon_model_type& __lexicon_source_target,
	      const lexicon_model_type& __lexicon_target_source)
    : LexiconBase(__lexicon_source_target, __lexicon_target_source) {}

  void assign_source(const utils::piece& phrase)
  {
    source.clear();
    cicada::TreeRule(phrase).frontier(std::back_inserter(source));
  }

  void assign_target(const utils::piece& phrase)
  {
    target.clear();
    cicada::TreeRule(phrase).frontier(std::back_inserter(target));
  }
};

#endif

