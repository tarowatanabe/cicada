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
#include <cicada/alignment.hpp>
#include <cicada/tree_rule.hpp>

struct RootCount
{
  typedef std::string label_type;
  typedef double count_type;
  typedef std::vector<count_type, std::allocator<count_type> > counts_type;
  
  label_type  label;
  counts_type counts;
  
  double observed;
  
  RootCount() : label(), counts() {}
  RootCount(const label_type& __label) : label(__label), counts() {}
  
  void clear()
  {
    label.clear();
    counts.clear();
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
  typedef std::string alignment_set_type;
  typedef double count_type;
  typedef std::vector<count_type, std::allocator<count_type> > counts_type;

  phrase_type source;
  phrase_type target;
  alignment_set_type alignments;
  counts_type counts;
  counts_type counts_source;
  counts_type counts_target;
  double observed_source;
  double observed_target;
  
  PhrasePair()
    : source(), target(), counts(), counts_source(), counts_target() {}

  void clear()
  {
    source.clear();
    target.clear();
    alignments.clear();
    counts.clear();
    counts_source.clear();
    counts_target.clear();
    observed_source = 0;
    observed_target = 0;
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
			  (RootCount::count_type, observed)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  PhrasePair,
			  (PhrasePair::phrase_type, source)
			  (PhrasePair::phrase_type, target)
			  (PhrasePair::alignment_set_type, alignments)
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
      root_count %= label >> "|||" >> counts >> "|||" >> qi::double_;
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
  
  typedef phrase_pair_type::phrase_type        phrase_type;
  typedef phrase_pair_type::alignment_set_type alignment_set_type;
  typedef phrase_pair_type::counts_type        counts_type;
  
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
      aligns %= qi::lexeme[+(standard::char_ - (standard::space >> "|||" >> standard::space))];
      counts %= +qi::double_;
      
      phrase_pair %= (phrase
		      >> "|||" >> phrase
		      >> "|||" >> aligns
		      >> "|||" >> counts
		      >> "|||" >> counts
		      >> "|||" >> counts
		      >> "|||" >> qi::double_ >> qi::double_);
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, alignment_set_type(), boost::spirit::standard::space_type> aligns;
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
  
  typedef phrase_pair_type::phrase_type        phrase_type;
  typedef phrase_pair_type::alignment_set_type alignment_set_type;
  typedef phrase_pair_type::counts_type        counts_type;
  
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
		      << " ||| " << standard::string
		      << " ||| " << counts
		      << " ||| " << counts
		      << " ||| " << counts
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

struct ExtractAlignment
{
  typedef std::pair<int, int> point_type;
  typedef cicada::Alignment alignment_type;
  typedef std::vector<alignment_type, std::allocator<alignment_type> > alignment_set_type;
  
  template <typename Iterator>
  struct align_parser : boost::spirit::qi::grammar<Iterator, alignment_set_type(), boost::spirit::standard::space_type>
  {
    align_parser() : align_parser::base_type(alignments)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      
      point %= qi::int_ >> '-' >> qi::int_;
      aligns %= '(' >> *point >> ')';
      alignments %= *aligns;
    }

    boost::spirit::qi::rule<Iterator, point_type(), boost::spirit::standard::space_type> point;
    boost::spirit::qi::rule<Iterator, alignment_type(), boost::spirit::standard::space_type> aligns;
    boost::spirit::qi::rule<Iterator, alignment_set_type(), boost::spirit::standard::space_type> alignments;
  };

  align_parser<std::string::const_iterator> parser;
  
  void operator()(const utils::piece& align, alignment_set_type& alignments)
  {
    namespace qi = boost::spirit::qi;
    namespace standard = boost::spirit::standard;
    
    std::string::const_iterator iter(align.begin());
    std::string::const_iterator end(align.end());
    
    alignments.clear();
    
    const bool result = qi::phrase_parse(iter, end, parser, standard::space, alignments);
    if (! result || iter != end)
      throw std::runtime_error("invalid point format? " + align);
  }
};

struct ExtractPhrase
{
  typedef cicada::Sentence  sentence_type;
  
  void operator()(const utils::piece& phrase, sentence_type& sentence)
  {
    sentence.assign(phrase);
  }
};

struct ExtractSCFG
{
  typedef cicada::Sentence  sentence_type;
  
  void operator()(const utils::piece& phrase, sentence_type& sentence)
  {
    sentence.assign(phrase);
    sentence.erase(sentence.begin());
  }
};

struct ExtractGHKM
{
  typedef cicada::Sentence  sentence_type;
  
  void operator()(const utils::piece& phrase, sentence_type& sentence)
  {
    sentence.clear();
    cicada::TreeRule(phrase).frontier(std::back_inserter(sentence));
  }
};


struct LexiconModel
{
  typedef cicada::Symbol word_type;
  typedef cicada::Vocab vocab_type;

  typedef boost::filesystem::path path_type;

  struct table_type
  {
    typedef utils::dense_hash_map<word_type, double, boost::hash<word_type> , std::equal_to<word_type> >::type __table_type;

    table_type() : table() { table.set_empty_key(word_type()); }

    __table_type table;
  };

  typedef utils::alloc_vector<table_type, std::allocator<table_type> > table_set_type;
  typedef std::vector<double, std::allocator<double> > max_set_type;

  LexiconModel() : smooth(1e-7), tables(), maximum() {}
  LexiconModel(const path_type& path) : smooth(), tables(), maximum() { open(path); }

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
    maximum.clear();
    
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
      
      const word_type source(boost::fusion::get<1>(parsed));
      const word_type target(boost::fusion::get<0>(parsed));
      
      tables[source.id()].table[target] = boost::fusion::get<2>(parsed);
      
      if (target.id() >= maximum.size())
	maximum.resize(target.id() + 1, 0.0);
      
      maximum[target.id()] = std::max(maximum[target.id()], boost::fusion::get<2>(parsed));
      
      smooth = std::min(smooth, boost::fusion::get<2>(parsed));
    }

    max_set_type(maximum).swap(maximum);
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
  
  double operator()(const word_type& word) const
  {
    return (word.id() < maximum.size() ?  std::max(smooth, maximum[word.id()]) : smooth);
  }

  double smooth;
  table_set_type tables;
  max_set_type   maximum;
};

struct Lexicon
{
  typedef cicada::Sentence  sentence_type;
  
  typedef ExtractAlignment::alignment_type     alignment_type;
  typedef ExtractAlignment::alignment_set_type alignment_set_type;

  typedef LexiconModel lexicon_model_type;

  typedef std::vector<int, std::allocator<int> > align_set_type;
  typedef std::vector<align_set_type, std::allocator<align_set_type> > align_map_type;

  Lexicon(const lexicon_model_type& __lexicon_source_target,
	  const lexicon_model_type& __lexicon_target_source)
    : lexicon_source_target(__lexicon_source_target),
      lexicon_target_source(__lexicon_target_source) {}

  const lexicon_model_type& lexicon_source_target;
  const lexicon_model_type& lexicon_target_source;

  align_map_type aligns_source_impl;
  align_map_type aligns_target_impl;
  
  std::pair<double, double> lexicon(const sentence_type& source,
				    const sentence_type& target,
				    const alignment_set_type& alignments) const
  {
    double lex_source_target = - std::numeric_limits<double>::infinity();
    double lex_target_source = - std::numeric_limits<double>::infinity();
    
    alignment_set_type::const_iterator aiter_end = alignments.end();
    for (alignment_set_type::const_iterator aiter = alignments.begin(); aiter != aiter_end; ++ aiter) {
      const std::pair<double, double> lex = lexicon(source, target, *aiter);
      
      lex_source_target = std::max(lex_source_target, lex.first);
      lex_target_source = std::max(lex_target_source, lex.second);
    }
    
    return std::make_pair(lex_source_target == - std::numeric_limits<double>::infinity() ? 0.0 : lex_source_target,
			  lex_target_source == - std::numeric_limits<double>::infinity() ? 0.0 : lex_target_source);
  }
  
  std::pair<double, double> lexicon(const sentence_type& source,
				    const sentence_type& target,
				    const alignment_type& alignment) const
  {
    const size_t source_size = source.size();
    const size_t target_size = target.size();

    double lex_source_target = 0.0;
    double lex_target_source = 0.0;

    align_map_type& aligns_source = const_cast<align_map_type&>(aligns_source_impl);
    align_map_type& aligns_target = const_cast<align_map_type&>(aligns_target_impl);

    aligns_source.clear();
    aligns_target.clear();

    aligns_source.resize(source_size);
    aligns_target.resize(target_size);
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      if (aiter->source >= static_cast<int>(source_size) || aiter->target >= static_cast<int>(target_size)) {
	std::ostringstream os;
	os << "invlaid alignment:"
	   << " source: " << source 
	   << " target: " << target
	   << " alignment: " << alignment;
	throw std::runtime_error(os.str());
      }
      
      aligns_source[aiter->source].push_back(aiter->target);
      aligns_target[aiter->target].push_back(aiter->source);
    }

    for (size_t trg = 0; trg != target_size; ++ trg) 
      if (target[trg].is_terminal()) {
	if (aligns_target[trg].empty())
	  lex_source_target += utils::mathop::log(lexicon_source_target(lexicon_model_type::vocab_type::EPSILON, target[trg]));
	else {
	  double score = 0.0;
	  align_set_type::const_iterator aiter_end = aligns_target[trg].end();
	  for (align_set_type::const_iterator aiter = aligns_target[trg].begin(); aiter != aiter_end; ++ aiter)
	    score += lexicon_source_target(source[*aiter], target[trg]);

	  lex_source_target += utils::mathop::log(score / aligns_target[trg].size());
	}
      }
    
    for (size_t src = 0; src != source_size; ++ src) 
      if (source[src].is_terminal()) {
	if (aligns_source[src].empty())
	  lex_target_source += utils::mathop::log(lexicon_target_source(lexicon_model_type::vocab_type::EPSILON, source[src]));
	else {
	  double score = 0.0;
	  align_set_type::const_iterator aiter_end = aligns_source[src].end();
	  for (align_set_type::const_iterator aiter = aligns_source[src].begin(); aiter != aiter_end; ++ aiter)
	    score += lexicon_target_source(target[*aiter], source[src]);
	  
	  lex_target_source += utils::mathop::log(score / aligns_source[src].size());
	}
      }
    
    return std::make_pair(lex_source_target, lex_target_source);
  }
  
  std::pair<double, double> noisy_or(const sentence_type& source,
				     const sentence_type& target) const
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
  
  std::pair<double, double> model1(const sentence_type& source,
				   const sentence_type& target) const
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
  
  std::pair<double, double> insertion_deletion(const sentence_type& source,
					       const sentence_type& target,
					       const double threshold_insertion=0.5,
					       const double threshold_deletion=0.5) const
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
	
	score_source_target -= (score < lexicon_source_target(target[trg]) * threshold_insertion);
      }
    
    double score_target_source = 0.0;
    
    for (size_t src = 0; src != source_size; ++ src)
      if (source[src].is_terminal()) {
	
	double score = 0.0;
	for (size_t trg = 0; trg != target_size; ++ trg)
	  if (target[trg].is_terminal())
	    score = std::max(score, lexicon_target_source(target[trg], source[src]));
	
	score_target_source -= (score < lexicon_target_source(source[src]) * threshold_deletion);
      }
    
    return std::make_pair(score_source_target, score_target_source);
  }
};

struct Unaligned
{
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef cicada::Sentence sentence_type;

  typedef ExtractAlignment::alignment_type     alignment_type;
  typedef ExtractAlignment::alignment_set_type alignment_set_type;
  
  typedef std::vector<int, std::allocator<int> > align_set_type;
  
  struct scores_type
  {
    size_type aligned_source;
    size_type aligned_target;
    size_type unaligned_source;
    size_type unaligned_target;
    
    scores_type() : aligned_source(0), aligned_target(0), unaligned_source(0), unaligned_target(0) {}
    scores_type(const size_type& __aligned_source,
		const size_type& __aligned_target,
		const size_type& __unaligned_source,
		const size_type& __unaligned_target)
      : aligned_source(__aligned_source),
        aligned_target(__aligned_target),
        unaligned_source(__unaligned_source),
        unaligned_target(__unaligned_target) {}
  };
  
  scores_type operator()(const sentence_type& source,
			 const sentence_type& target,
			 const alignment_set_type& alignments) const
  {
    scores_type scores(0, 0, size_type(-1), size_type(-1));
    
    alignment_set_type::const_iterator aiter_end = alignments.end();
    for (alignment_set_type::const_iterator aiter = alignments.begin(); aiter != aiter_end; ++ aiter) {
      const scores_type result = operator()(source, target, *aiter);
      
      scores.aligned_source = utils::bithack::max(scores.aligned_source, result.aligned_source);
      scores.aligned_target = utils::bithack::max(scores.aligned_target, result.aligned_target);
      
      scores.unaligned_source = utils::bithack::min(scores.unaligned_source, result.unaligned_source);
      scores.unaligned_target = utils::bithack::min(scores.unaligned_target, result.unaligned_target);
    }
    
    return scores_type(scores.aligned_source,
		       scores.aligned_target,
		       utils::bithack::branch(scores.unaligned_source == size_type(-1), size_type(0), scores.unaligned_source),
		       utils::bithack::branch(scores.unaligned_target == size_type(-1), size_type(0), scores.unaligned_target));
  }
  
  align_set_type aligns_source_impl;
  align_set_type aligns_target_impl;
  
  scores_type operator()(const sentence_type& source,
			 const sentence_type& target,
			 const alignment_type& alignment) const
  {
    const size_t source_size = source.size();
    const size_t target_size = target.size();
    
    align_set_type& aligns_source = const_cast<align_set_type&>(aligns_source_impl);
    align_set_type& aligns_target = const_cast<align_set_type&>(aligns_target_impl);
    
    aligns_source.clear();
    aligns_target.clear();
    
    aligns_source.resize(source_size);
    aligns_target.resize(target_size);
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      if (aiter->source >= static_cast<int>(source_size) || aiter->target >= static_cast<int>(target_size)) {
	std::ostringstream os;
	os << "invlaid alignment:"
	   << " source: " << source 
	   << " target: " << target
	   << " alignment: " << alignment;
	throw std::runtime_error(os.str());
      }
      
      ++ aligns_source[aiter->source];
      ++ aligns_target[aiter->target];
    }
    
    size_type terminals_source = 0;
    size_type terminals_target = 0;
    size_type unaligned_source = 0;
    size_type unaligned_target = 0;
    
    for (size_t src = 0; src != source_size; ++ src) 
      if (source[src].is_terminal()) {
	unaligned_source += (aligns_source[src] == 0);
	++ terminals_source;
      }
    
    for (size_t trg = 0; trg != target_size; ++ trg) 
      if (target[trg].is_terminal()) {
	unaligned_target += (aligns_target[trg] == 0);
	++ terminals_target;
      }
    
    return scores_type(terminals_source - unaligned_source, 
		       terminals_target - unaligned_target, 
		       unaligned_source,
		       unaligned_target);
  }
  
};

struct Cross
{
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef cicada::Sentence sentence_type;

  typedef ExtractAlignment::alignment_type     alignment_type;
  typedef ExtractAlignment::alignment_set_type alignment_set_type;

  typedef std::vector<size_type, std::allocator<size_type> > position_type;

  size_type operator()(const sentence_type& source,
		       const sentence_type& target,
		       const alignment_set_type& alignments)
  {
    size_type crossed(size_type(-1));
    
    alignment_set_type::const_iterator aiter_end = alignments.end();
    for (alignment_set_type::const_iterator aiter = alignments.begin(); aiter != aiter_end; ++ aiter) {
      const size_type result = operator()(source, target, *aiter);
      
      crossed = utils::bithack::min(crossed, result);
    }
    
    return utils::bithack::branch(crossed == size_type(-1), size_type(0), crossed);
  }
  
  size_type operator()(const sentence_type& source,
		       const sentence_type& target,
		       const alignment_type& alignment)
  {
    if (alignment.empty()) return 0;

    size_type crossed = 0;
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end - 1; ++ aiter)
      for (alignment_type::const_iterator niter = aiter + 1; niter != aiter_end; ++ niter)
	crossed += ((aiter->source < niter->source && niter->target < aiter->target)
		    || (niter->source < aiter->source && aiter->target < niter->target));
    
    return crossed;
  }
  
  size_type operator()(const sentence_type& source,
		       const sentence_type& target)
  {
    //
    // compute non-terminal alignment
    //
    
    // first, enumerate source-position...
    position_source.clear();
    
    int pos = 1;
    sentence_type::const_iterator siter_begin = source.begin();
    sentence_type::const_iterator siter_end   = source.end();
    for (sentence_type::const_iterator siter = siter_begin; siter != siter_end; ++ siter)
      if (siter->is_non_terminal()) {
	const int non_terminal_pos = siter->non_terminal_index();
	const int non_terminal_index = utils::bithack::branch(non_terminal_pos <= 0, pos, non_terminal_pos);
	
	if (non_terminal_index >= static_cast<int>(position_source.size()))
	  position_source.resize(non_terminal_index + 1);
	
	position_source[non_terminal_index] = siter - siter_begin;
	
	++ pos;
      }
    
    if (pos == 1) return 0;
    
    // then, enumerate target-position...
    alignment.clear();

    pos = 1;
    sentence_type::const_iterator titer_begin = target.begin();
    sentence_type::const_iterator titer_end   = target.end();
    for (sentence_type::const_iterator titer = titer_begin; titer != titer_end; ++ titer)
      if (titer->is_non_terminal()) {
	const int non_terminal_pos = titer->non_terminal_index();
	const int non_terminal_index = utils::bithack::branch(non_terminal_pos <= 0, pos, non_terminal_pos);
	
	alignment.push_back(alignment_type::value_type(position_source[non_terminal_index], titer - titer_begin));
	
	++ pos;
      }
    
    //
    // compute crossing...
    //
    
    size_type crossed = 0;
    
    // check whether niter crossed agains aiter...
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end - 1; ++ aiter)
      for (alignment_type::const_iterator niter = aiter + 1; niter != aiter_end; ++ niter)
	crossed += ((aiter->source < niter->source && niter->target < aiter->target)
		    || (niter->source < aiter->source && aiter->target < niter->target));
    
    return crossed;
  }
  
  alignment_type alignment;
  position_type  position_source;
};


#endif

