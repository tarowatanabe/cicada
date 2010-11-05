#ifndef __CICADA__EXTRACT_SCORE_IMPL__HPP__
#define __CICADA__EXTRACT_SCORE_IMPL__HPP__ 1

#define BOOST_SPIRIT_THREADSAFE
#define PHOENIX_THREADSAFE

#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/karma.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>

#include <boost/fusion/tuple.hpp>
#include <boost/fusion/adapted.hpp>

#include <boost/tuple/tuple.hpp>
#include <boost/filesystem.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/array.hpp>

#include <string>
#include <vector>
#include <set>
#include <algorithm>
#include <utility>
#include <queue>
#include <deque>
#include <sstream>
#include <iostream>

#include <cicada/symbol.hpp>
#include <cicada/vocab.hpp>
#include <cicada/sentence.hpp>
#include <cicada/tree_rule.hpp>

#include <utils/tempfile.hpp>
#include <utils/hashmurmur.hpp>
#include <utils/lockfree_list_queue.hpp>
#include <utils/compress_stream.hpp>
#include <utils/array_power2.hpp>
#include <utils/repository.hpp>
#include <utils/map_file.hpp>
#include <utils/alloc_vector.hpp>
#include <utils/sgi_hash_set.hpp>
#include <utils/malloc_stats.hpp>
#include <utils/lexical_cast.hpp>

#include <succinct_db/succinct_trie_db.hpp>

#include <google/dense_hash_map>

struct RootCount
{
  typedef std::string label_type;
  typedef std::vector<double, std::allocator<double> > counts_type;
  
  label_type  label;
  counts_type counts;

  RootCount(const label_type& __label) : label(__label), counts() {}
  RootCount() : label(), counts() {}

  void clear()
  {
    label.clear();
    counts.clear();
  }

  void swap(RootCount& x)
  {
    label.swap(x.label);
    counts.swap(x.counts);
  }

  
  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    counts.resize(last - first, 0.0);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
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
  typedef std::vector<double, std::allocator<double> > counts_type;
  typedef std::pair<int, int> point_type;
  typedef std::vector<point_type, std::allocator<point_type> > alignment_type;

  phrase_type    source;
  phrase_type    target;
  alignment_type alignment;
  counts_type    counts;

  PhrasePair() : source(), target(), alignment(), counts() {}

  void clear()
  {
    source.clear();
    target.clear();
    alignment.clear();
    counts.clear();
  }

  void swap(PhrasePair& x)
  {
    source.swap(x.source);
    target.swap(x.target);
    alignment.swap(x.alignment);
    counts.swap(x.counts);
  }
  
  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    counts.resize(last - first, 0.0);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
  }

  friend
  size_t  hash_value(PhrasePair const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;
    
    return hasher_type()(x.source.begin(), x.source.end(),
			 hasher_type()(x.target.begin(), x.target.end(),
				       hasher_type()(x.alignment.begin(), x.alignment.end(), 0)));
  }

  friend
  bool operator==(const PhrasePair& x, const PhrasePair& y) 
  {
    return x.source == y.source && x.target == y.target && x.alignment == y.alignment;
  }
  
  friend
  bool operator!=(const PhrasePair& x, const PhrasePair& y) 
  {
    return x.source != y.source || x.target != y.target || x.alignment != y.alignment;
  }
  
  friend
  bool operator<(const PhrasePair& x, const PhrasePair& y)
  {
    return (x.source < y.source
	    || (!(y.source < x.source)
		&& (x.target < y.target
		    || (!(y.target < x.target)
			&& x.alignment < y.alignment))));
  }

  friend
  bool operator>(const PhrasePair& x, const PhrasePair& y)
  {
    return y < x;
  }
};

struct PhrasePairModified
{
  typedef PhrasePair phrase_pair_type;
  typedef phrase_pair_type::phrase_type phrase_type;
  typedef phrase_pair_type::counts_type counts_type;

  phrase_type    source;
  phrase_type    target;
  counts_type    counts;

  PhrasePairModified() : source(), target(), counts() {}

  void clear()
  {
    source.clear();
    target.clear();
    counts.clear();
  }

  void swap(PhrasePairModified& x)
  {
    source.swap(x.source);
    target.swap(x.target);
    counts.swap(x.counts);
  }
  
  template <typename Iterator>
  void increment(Iterator first, Iterator last)
  {
    counts.resize(last - first);
    std::transform(first, last, counts.begin(), counts.begin(), std::plus<double>());
  }
  
  friend
  size_t  hash_value(PhrasePairModified const& x)
  {
    typedef utils::hashmurmur<size_t> hasher_type;
   
    return hasher_type()(x.source.begin(), x.source.end(), hasher_type()(x.target.begin(), x.target.end(),0));
  }


  friend
  bool operator==(const PhrasePairModified& x, const PhrasePairModified& y) 
  {
    return x.source == y.source && x.target == y.target;
  }
  
  friend
  bool operator!=(const PhrasePairModified& x, const PhrasePairModified& y) 
  {
    return x.source != y.source || x.target != y.target;
  }
  
  friend
  bool operator<(const PhrasePairModified& x, const PhrasePairModified& y)
  {
    return (x.source < y.source || (!(y.source < x.source) && x.target < y.target));
  }

  friend
  bool operator>(const PhrasePairModified& x, const PhrasePairModified& y)
  {
    return y < x;
  }
  
};

namespace std
{
  inline
  void swap(RootCount& x, RootCount& y)
  {
    x.swap(y);
  }
  
  inline
  void swap(PhrasePair& x, PhrasePair& y)
  {
    x.swap(y);
  }

  inline
  void swap(PhrasePairModified& x, PhrasePairModified& y)
  {
    x.swap(y);
  }
};

BOOST_FUSION_ADAPT_STRUCT(
			  RootCount,
			  (RootCount::label_type, label)
			  (RootCount::counts_type, counts)
			  )

BOOST_FUSION_ADAPT_STRUCT(
			  PhrasePair,
			  (PhrasePair::phrase_type, source)
			  (PhrasePair::phrase_type, target)
			  (PhrasePair::alignment_type, alignment)
			  (PhrasePair::counts_type, counts)
			  )
BOOST_FUSION_ADAPT_STRUCT(
			  PhrasePairModified,
			  (PhrasePairModified::phrase_type, source)
			  (PhrasePairModified::phrase_type, target)
			  (PhrasePairModified::counts_type, counts)
			  )


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
    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();

    std::string label;

    const bool result = boost::spirit::qi::phrase_parse(iter, end,
							boost::spirit::qi::lexeme[+(boost::spirit::standard::char_ - boost::spirit::standard::space)],
							boost::spirit::standard::space,
							label);
    if (! result)
      throw std::runtime_error("no label?");
    
    return label;
  }
};

struct ExtractRootTree
{
  // extract the first word...
  std::string operator()(const std::string& phrase) const
  {
    std::string::const_iterator iter = phrase.begin();
    std::string::const_iterator end = phrase.end();
    
    std::string label;
    
    const bool result = boost::spirit::qi::phrase_parse(iter, end,
							boost::spirit::qi::lexeme[+((boost::spirit::standard::char_ - boost::spirit::standard::space - '(') | "\\(")],
							boost::spirit::standard::space,
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
    typedef boost::fusion::tuple<std::string, std::string, double> parsed_type;

    boost::spirit::qi::rule<std::string::const_iterator, std::string(), boost::spirit::standard::space_type> word;
    boost::spirit::qi::rule<std::string::const_iterator, parsed_type(), boost::spirit::standard::space_type> lexicon;
    
    word %= boost::spirit::qi::lexeme[+(boost::spirit::standard::char_ - boost::spirit::standard::space)];
    lexicon %= word >> word >> boost::spirit::qi::double_;
    
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
      
      const bool result = boost::spirit::qi::phrase_parse(iter, end, lexicon, boost::spirit::standard::space, parsed);
      if (! result || iter != end) continue;
      
      tables[word_type(boost::fusion::get<1>(parsed)).id()].table[boost::fusion::get<0>(parsed)] = boost::fusion::get<2>(parsed);
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

struct LexiconPhrase
{
  typedef PhrasePair::point_type     point_type;
  typedef PhrasePair::alignment_type alignment_type;
  
  typedef cicada::Sentence sentence_type;

  typedef std::vector<int, std::allocator<int> > align_set_type;
  typedef std::vector<align_set_type, std::allocator<align_set_type> > align_map_type;
  
  typedef LexiconModel lexicon_model_type;

  LexiconPhrase(const lexicon_model_type& __lexicon_source_target,
		const lexicon_model_type& __lexicon_target_source)
    : lexicon_source_target(__lexicon_source_target),
      lexicon_target_source(__lexicon_target_source) {}
  
  const lexicon_model_type& lexicon_source_target;
  const lexicon_model_type& lexicon_target_source;

  align_map_type aligns_source_impl;
  align_map_type aligns_target_impl;
  
  std::pair<double, double> operator()(const std::string& source_str, const std::string& target_str, const alignment_type& alignment) const
  {
    const sentence_type source(source_str);
    const sentence_type target(target_str);

    const size_t source_size = source.size();
    const size_t target_size = target.size();
    
    double lex_source_target = 1.0;
    double lex_target_source = 1.0;

    align_map_type& aligns_source = const_cast<align_map_type&>(aligns_source_impl);
    align_map_type& aligns_target = const_cast<align_map_type&>(aligns_target_impl);
    
    aligns_source.clear();
    aligns_target.clear();
    
    aligns_source.resize(source_size);
    aligns_target.resize(target_size);
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      aligns_source[aiter->first].push_back(aiter->second);
      aligns_target[aiter->second].push_back(aiter->first);
    }
    
    for (size_t trg = 0; trg != target_size; ++ trg) {
      if (aligns_target[trg].empty())
	lex_source_target *= lexicon_source_target(lexicon_model_type::vocab_type::NONE, target[trg]);
      else {
	double score = 0.0;
	align_set_type::const_iterator aiter_end = aligns_target[trg].end();
	for (align_set_type::const_iterator aiter = aligns_target[trg].begin(); aiter != aiter_end; ++ aiter)
	  score += lexicon_source_target(source[*aiter], target[trg]);
	
	lex_source_target *= score / aligns_target[trg].size();
      }
    }
    
    for (size_t src = 0; src != source_size; ++ src) {
      if (aligns_source[src].empty())
	lex_target_source *= lexicon_target_source(lexicon_model_type::vocab_type::NONE, source[src]);
      else {
	double score = 0.0;
	align_set_type::const_iterator aiter_end = aligns_source[src].end();
	for (align_set_type::const_iterator aiter = aligns_source[src].begin(); aiter != aiter_end; ++ aiter)
	  score += lexicon_target_source(target[*aiter], source[src]);
	
	lex_target_source *= score / aligns_source[src].size();
      }
    }
    
    return std::make_pair(lex_source_target, lex_target_source);
  }
};

struct LexiconSCFG
{
  typedef PhrasePair::point_type     point_type;
  typedef PhrasePair::alignment_type alignment_type;
  
  typedef cicada::Sentence sentence_type;

  typedef std::vector<int, std::allocator<int> > align_set_type;
  typedef std::vector<align_set_type, std::allocator<align_set_type> > align_map_type;
  
  typedef LexiconModel lexicon_model_type;

  LexiconSCFG(const lexicon_model_type& __lexicon_source_target,
	      const lexicon_model_type& __lexicon_target_source)
    : lexicon_source_target(__lexicon_source_target),
      lexicon_target_source(__lexicon_target_source) {}
  
  const lexicon_model_type& lexicon_source_target;
  const lexicon_model_type& lexicon_target_source;
  
  align_map_type aligns_source_impl;
  align_map_type aligns_target_impl;
  
  std::pair<double, double> operator()(const std::string& source_str, const std::string& target_str, const alignment_type& alignment) const
  {
    // major difference to LexiconPhrase is the addition of root-symbol (but we will keep the same alignment structure...)
    
    const sentence_type source(source_str);
    const sentence_type target(target_str);

    const size_t source_size = source.size() - 1;
    const size_t target_size = target.size() - 1;
    
    double lex_source_target = 1.0;
    double lex_target_source = 1.0;

    align_map_type& aligns_source = const_cast<align_map_type&>(aligns_source_impl);
    align_map_type& aligns_target = const_cast<align_map_type&>(aligns_target_impl);
    
    aligns_source.clear();
    aligns_target.clear();
    
    aligns_source.resize(source_size);
    aligns_target.resize(target_size);
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      aligns_source[aiter->first].push_back(aiter->second);
      aligns_target[aiter->second].push_back(aiter->first);
    }
    
    for (size_t trg = 0; trg != target_size; ++ trg) 
      if (target[trg].is_terminal()) {
	if (aligns_target[trg].empty())
	  lex_source_target *= lexicon_source_target(lexicon_model_type::vocab_type::NONE, target[trg + 1]);
	else {
	  double score = 0.0;
	  align_set_type::const_iterator aiter_end = aligns_target[trg].end();
	  for (align_set_type::const_iterator aiter = aligns_target[trg].begin(); aiter != aiter_end; ++ aiter)
	    score += lexicon_source_target(source[*aiter + 1], target[trg + 1]);
	  
	  lex_source_target *= score / aligns_target[trg].size();
	}
      }
    
    for (size_t src = 0; src != source_size; ++ src) 
      if (source[src].is_terminal()) {
	if (aligns_source[src].empty())
	  lex_target_source *= lexicon_target_source(lexicon_model_type::vocab_type::NONE, source[src + 1]);
	else {
	  double score = 0.0;
	  align_set_type::const_iterator aiter_end = aligns_source[src].end();
	  for (align_set_type::const_iterator aiter = aligns_source[src].begin(); aiter != aiter_end; ++ aiter)
	    score += lexicon_target_source(target[*aiter + 1], source[src + 1]);
	  
	  lex_target_source *= score / aligns_source[src].size();
	}
      }
    
    return std::make_pair(lex_source_target, lex_target_source);
  }
};

struct LexiconTree
{
  typedef PhrasePair::point_type     point_type;
  typedef PhrasePair::alignment_type alignment_type;
  
  typedef cicada::Sentence sentence_type;
  typedef cicada::TreeRule rule_type;

  typedef std::vector<int, std::allocator<int> > align_set_type;
  typedef std::vector<align_set_type, std::allocator<align_set_type> > align_map_type;
  
  typedef LexiconModel lexicon_model_type;

  LexiconTree(const lexicon_model_type& __lexicon_source_target,
	      const lexicon_model_type& __lexicon_target_source)
    : lexicon_source_target(__lexicon_source_target),
      lexicon_target_source(__lexicon_target_source) {}

  const lexicon_model_type& lexicon_source_target;
  const lexicon_model_type& lexicon_target_source;
  
  align_map_type aligns_source_impl;
  align_map_type aligns_target_impl  ;
  
  
  
  std::pair<double, double> operator()(const std::string& source_str, const std::string& target_str, const alignment_type& alignment) const
  {
    // we will parse only the frontiers from source-str and target-str!
    
    sentence_type source;
    sentence_type target;

    rule_type(source_str).frontier(std::back_inserter(source));
    rule_type(target_str).frontier(std::back_inserter(target));
    
    const size_t source_size = source.size();
    const size_t target_size = target.size();
    
    double lex_source_target = 1.0;
    double lex_target_source = 1.0;
    
    align_map_type& aligns_source = const_cast<align_map_type&>(aligns_source_impl);
    align_map_type& aligns_target = const_cast<align_map_type&>(aligns_target_impl);

    aligns_source.clear();
    aligns_target.clear();
    
    aligns_source.resize(source_size);
    aligns_target.resize(target_size);
    
    alignment_type::const_iterator aiter_end = alignment.end();
    for (alignment_type::const_iterator aiter = alignment.begin(); aiter != aiter_end; ++ aiter) {
      aligns_source[aiter->first].push_back(aiter->second);
      aligns_target[aiter->second].push_back(aiter->first);
    }
    
    for (size_t trg = 0; trg != target_size; ++ trg) 
      if (target[trg].is_terminal()) {
	if (aligns_target[trg].empty())
	  lex_source_target *= lexicon_source_target(lexicon_model_type::vocab_type::NONE, target[trg]);
	else {
	  double score = 0.0;
	  align_set_type::const_iterator aiter_end = aligns_target[trg].end();
	  for (align_set_type::const_iterator aiter = aligns_target[trg].begin(); aiter != aiter_end; ++ aiter)
	    score += lexicon_source_target(source[*aiter], target[trg]);
	  
	  lex_source_target *= score / aligns_target[trg].size();
	}
      }
    
    for (size_t src = 0; src != source_size; ++ src) 
      if (source[src].is_terminal()) {
	if (aligns_source[src].empty())
	  lex_target_source *= lexicon_target_source(lexicon_model_type::vocab_type::NONE, source[src]);
	else {
	  double score = 0.0;
	  align_set_type::const_iterator aiter_end = aligns_source[src].end();
	  for (align_set_type::const_iterator aiter = aligns_source[src].begin(); aiter != aiter_end; ++ aiter)
	    score += lexicon_target_source(target[*aiter], source[src]);
	  
	  lex_target_source *= score / aligns_source[src].size();
	}
      }
    
    return std::make_pair(lex_source_target, lex_target_source);
  }
};


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
      namespace phoenix = boost::phoenix;
      
      using qi::phrase_parse;
      using qi::lexeme;
      using qi::hold;
      using standard::char_;
      using qi::double_;
      using standard::space;
      
      label %= lexeme[+(char_ - (space >> "|||" >> space))];
      counts %= +double_;
      root_count %= label >> "|||" >> counts;
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


struct RootCountGenerator
{
  typedef RootCount root_count_type;
  
  typedef root_count_type::label_type  label_type;
  typedef root_count_type::counts_type counts_type;
 
  RootCountGenerator() : grammar() {}
  RootCountGenerator(const RootCountGenerator& x) : grammar() {}
  
  template <typename Iterator>
  struct root_count_generator : boost::spirit::karma::grammar<Iterator, root_count_type()>
  {
    root_count_generator() : root_count_generator::base_type(root_count)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::repeat;
      using standard::char_;
      using karma::double_;
      using standard::space;
      
      label  %= +char_;
      counts %= double20 % ' ';
      root_count %= label << " ||| " << counts;
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 20;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double20;
    
    boost::spirit::karma::rule<Iterator, std::string()> label;
    boost::spirit::karma::rule<Iterator, counts_type()> counts;
    boost::spirit::karma::rule<Iterator, root_count_type()> root_count;
  };

  typedef std::ostream_iterator<char> iterator_type;
  
  std::ostream& operator()(std::ostream& os, const root_count_type& root_count)
  {
    iterator_type iter(os);
    
    if (! boost::spirit::karma::generate(iter, grammar, root_count))
      throw std::runtime_error("failed generation!");
    
    return os;
  }
  
  root_count_generator<iterator_type> grammar;
};

struct PhrasePairParser
{
  typedef PhrasePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
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
      namespace phoenix = boost::phoenix;
      
      using qi::lexeme;
      using qi::hold;
      using qi::repeat;
      using standard::char_;
      using qi::double_;
      using qi::int_;
      using standard::space;
      
      phrase %= lexeme[+(char_ - (space >> "|||" >> space))];
      alignment %= *(int_ >> '-' >> int_);
      counts %= +double_;
      phrase_pair %= phrase >> "|||" >> phrase >> "|||" >> alignment >> "|||" >> counts;
    }
    
    boost::spirit::qi::rule<Iterator, std::string(), boost::spirit::standard::space_type> phrase;
    boost::spirit::qi::rule<Iterator, alignment_type(), boost::spirit::standard::space_type> alignment;
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
    
    const bool result =  boost::spirit::qi::phrase_parse(iter, end, grammar, boost::spirit::standard::space, phrase_pair);
    return result && iter == end;
  }
  
  phrase_pair_parser<std::string::const_iterator> grammar;
};


struct PhrasePairGenerator
{
  typedef PhrasePair phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::alignment_type alignment_type;
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
      namespace phoenix = boost::phoenix;
      
      using karma::repeat;
      using standard::char_;
      using karma::double_;
      using karma::int_;
      using standard::space;
      
      phrase %= +char_;
      alignment %= -((int_ << '-' << int_) % ' ');
      counts %= double20 % ' ';
      phrase_pair %= phrase << " ||| " << phrase << " ||| " << alignment << " ||| " << counts;
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 20;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double20;
    
    boost::spirit::karma::rule<Iterator, std::string()> phrase;
    boost::spirit::karma::rule<Iterator, alignment_type()> alignment;
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


struct PhrasePairModifiedParser
{
  typedef PhrasePairModified phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;
  
  PhrasePairModifiedParser() : grammar() {}
  PhrasePairModifiedParser(const PhrasePairModifiedParser& x) : grammar() {}
  
  template <typename Iterator>
  struct phrase_pair_parser : boost::spirit::qi::grammar<Iterator, phrase_pair_type(), boost::spirit::standard::space_type>
  {
    phrase_pair_parser() : phrase_pair_parser::base_type(phrase_pair)
    {
      namespace qi = boost::spirit::qi;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using qi::lexeme;
      using qi::hold;
      using standard::char_;
      using qi::double_;
      using standard::space;
      
      phrase %= lexeme[+(char_ - (space >> "|||" >> space))];
      counts %= +double_;
      phrase_pair %= phrase >> "|||" >> phrase >> "|||" >> counts;
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
    return result && iter == end;
  }
  
  phrase_pair_parser<std::string::const_iterator> grammar;
};

struct PhrasePairModifiedGenerator
{
  typedef PhrasePairModified phrase_pair_type;
  
  typedef phrase_pair_type::phrase_type    phrase_type;
  typedef phrase_pair_type::counts_type    counts_type;

  PhrasePairModifiedGenerator() : grammar() {}
  PhrasePairModifiedGenerator(const PhrasePairModifiedGenerator& x) : grammar() {}
  
  
  template <typename Iterator>
  struct phrase_pair_generator : boost::spirit::karma::grammar<Iterator, phrase_pair_type()>
  {
    phrase_pair_generator() : phrase_pair_generator::base_type(phrase_pair)
    {
      namespace karma = boost::spirit::karma;
      namespace standard = boost::spirit::standard;
      namespace phoenix = boost::phoenix;
      
      using karma::repeat;
      using standard::char_;
      using karma::double_;
      using standard::space;
      
      phrase %= +char_;
      counts %= double20 % ' ';
      phrase_pair %= phrase << " ||| " << phrase << " ||| " << " ||| " << counts;
    }

    struct real_precision : boost::spirit::karma::real_policies<double>
    {
      static unsigned int precision(double) 
      { 
        return 20;
      }
    };
    
    boost::spirit::karma::real_generator<double, real_precision> double20;
    
    boost::spirit::karma::rule<Iterator, std::string()> phrase;
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


// modify counts... simply map from source to target, and collect counts
struct PhrasePairModify
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;

  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;
  
  typedef RootCount          root_count_type;
  typedef PhrasePair         phrase_pair_type;
  typedef PhrasePairModified modified_type;
  
  typedef std::set<root_count_type, std::less<root_count_type>, std::allocator<root_count_type> > root_count_set_type;
  typedef std::vector<modified_type, std::allocator<modified_type> >  modified_set_type;
  
  typedef utils::lockfree_list_queue<modified_set_type, std::allocator<modified_set_type> > queue_type;
  
  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;

  
};

template <typename ExtractRoot>
struct PhrasePairModifyMapper
{
  typedef PhrasePairModify map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  
  typedef map_reduce_type::modified_type     modified_type;
  typedef map_reduce_type::modified_set_type modified_set_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef PhrasePairParser phrase_pair_parser_type;
  typedef ExtractRoot extract_root_type;
  
  hasher_type             hasher;
  phrase_pair_parser_type phrase_pair_parser;
  
  path_set_type paths;
  queue_ptr_set_type& queues;
  root_count_set_type& root_counts;
  const extract_root_type& extract_root;
  int debug;

  PhrasePairModifyMapper(const path_set_type& __paths,
			 queue_ptr_set_type& __queues,
			 root_count_set_type& __root_counts,
			 const extract_root_type& __extract_root,
			 const int __debug)
    : paths(__paths), queues(__queues), root_counts(__root_counts), extract_root(__extract_root), debug(__debug) {}

  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first.front() > y->first.front();
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first.front() > y->first.front();
    }
  };
  
  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    phrase_pair_type phrase_pair;
    
    while (counts.size() < 256 && std::getline(is, line)) {
      if (! phrase_pair_parser(line, phrase_pair)) continue;
      
      if (counts.empty() || counts.back().first.source != phrase_pair.source || counts.back().first.target != phrase_pair.target) {
	counts.resize(counts.size() + 1);
	counts.back().source.swap(phrase_pair.source);
	counts.back().target.swap(phrase_pair.target);
	counts.back().counts.swap(phrase_pair.counts);
      } else
	counts.back().increment(phrase_pair.counts.begin(), phrase_pair.counts.end());
    }
  }
  
  void operator()()
  {
    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
    
    typedef std::deque<modified_type, std::allocator<modified_type> > buffer_type;
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef boost::shared_ptr<buffer_stream_type> buffer_stream_ptr_type;
    typedef std::vector<buffer_stream_ptr_type, std::allocator<buffer_stream_ptr_type> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_ptr_type, pqueue_base_type, greater_buffer<buffer_stream_ptr_type> > pqueue_type;
    
    typedef std::vector<modified_set_type, std::allocator<modified_set_type> > modified_map_type;
    
    pqueue_type          pqueue;
    istream_ptr_set_type istreams;
    
    path_set_type::const_iterator piter_end = paths.end();
    for (path_set_type::const_iterator piter = paths.begin(); piter != piter_end; ++ piter) {
      if (! boost::filesystem::exists(*piter))
	throw std::runtime_error("no file? " + piter->file_string());

      istreams.push_back(istream_ptr_type(new istream_type(*piter, 1024 * 1024)));
      
      buffer_stream_ptr_type buffer_stream(new buffer_stream_type());
      buffer_stream->second = &(*istreams.back());
      
      read_phrase_pair(*istreams.back(), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    modified_map_type modified(queues.size());
    modified_type counts;

    int iter = 0;
    const int mask_full    = (1 << 12) - 1;
    const int mask_partial = (1 << 10) - 1;
    
    while (! pqueue.empty()) {
      buffer_stream_ptr_type buffer_stream(pqueue.top());
      pqueue.pop();
      
      if (counts != buffer_stream->first.front()) {
	if (! counts.counts.empty()) {
	  std::pair<root_count_set_type::iterator, bool> result_root = root_counts.insert(extract_root(counts.source));
	  if (! result_root.second)
	    const_cast<root_count_type&>(*result_root.first).increment(counts.counts.begin(), counts.counts.end());
	  
	  // swap source and target!
	  
	  counts.source.swap(counts.target);
	  const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
	  modified[shard].push_back(counts);
	}
	
	if ((iter & mask_full) == mask_full) {
	  for (int shard = 0; shard < queues.size(); ++ shard)
	    if (modified[shard].size() >= 256) {
	      queues[shard]->push_swap(modified[shard]);
	      modified[shard].clear();
	    }
	} else if ((iter & mask_partial == mask_partial)) {
	  bool found = false;
	  for (int shard = 0; shard < queues.size(); ++ shard)
	    if (modified[shard].size() >= 256) {
	      const bool no_wait = (modified[shard].size() < 1024 * 4);
	      if (queues[shard]->push_swap(modified[shard], no_wait)) {
		modified[shard].clear();
		found = true;
	      }
	    }
	  
	  if (! found)
	    boost::thread::yield();
	}
	
	++ iter;
	
	counts.swap(buffer_stream->first.front());
      } else
	counts.increment(buffer_stream->first.front().counts.begin(), buffer_stream->first.front().counts.end());
      
      buffer_stream->first.pop_front();
      
      if (buffer_stream->first.empty())
	read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    if (! counts.counts.empty()) {
      std::pair<root_count_set_type::iterator, bool> result_root = root_counts.insert(extract_root(counts.source));
      if (! result_root.second)
	const_cast<root_count_type&>(*result_root.first).increment(counts.counts.begin(), counts.counts.end());

      // swap source and target!
      counts.source.swap(counts.target);
      const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
      modified[shard].push_back(counts);
    }
    
    // send remaining...
    for (int shard = 0; shard < queues.size(); ++ shard)
      if (! modified[shard].empty()) {
	queues[shard]->push_swap(modified[shard]);
	modified[shard].clear();
      }
    
    // termination...
    
    while (1) {
      for (size_t shard = 0; shard != queues.size(); ++ shard) 
	if (queues[shard]) {
	  modified[shard].clear();
	  
	  if (queues[shard]->push_swap(modified[shard], true)) {
	    queues[shard].reset();
	    modified[shard].clear();
	  }
	}
      
      if (std::count(queues.begin(), queues.end(), queue_ptr_type()) == queues.size()) break;
      
      boost::thread::yield();
    }
  }
};

template <typename ExtractRoot>
struct PhrasePairModifyReducer
{
  typedef PhrasePairModify map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;
  
  typedef map_reduce_type::modified_type     modified_type;
  typedef map_reduce_type::modified_set_type modified_set_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;
  
  typedef PhrasePairModifiedParser    modified_parser_type;
  typedef PhrasePairModifiedGenerator modified_generator_type;

  typedef ExtractRoot extract_root_type;


#ifdef HAVE_TR1_UNORDERED_SET
  typedef std::tr1::unordered_set<modified_type, boost::hash<modified_type>, std::equal_to<modified_type>,
				  std::allocator<modified_type> > modified_unique_type;
#else
  typedef sgi::hash_set<modified_type, boost::hash<modified_type>, std::equal_to<modified_type>,
			std::allocator<modified_type> > modified_unique_type;
#endif
  
  modified_parser_type    parser;
  modified_generator_type generator;

  queue_type&    queue;
  path_set_type& paths;
  root_count_set_type& root_counts;
  const extract_root_type& extract_root;
  
  int            shard_size;
  double         max_malloc;
  int            debug;
  
  PhrasePairModifyReducer(queue_type&    __queue,
			  path_set_type& __paths,
			  root_count_set_type& __root_counts,
			  const extract_root_type& __extract_root,
			  const int      __shard_size,
			  const double   __max_malloc,
			  const int      __debug)
    : queue(__queue),
      paths(__paths),
      root_counts(__root_counts),
      extract_root(__extract_root),
      shard_size(__shard_size),
      max_malloc(__max_malloc),
      debug(__debug) {}

  struct less_file_size
  {
    bool operator()(const path_type& x, const path_type& y) const
    {
      return boost::filesystem::file_size(x) < boost::filesystem::file_size(y);
    }
  };


  // merge counts from two streams into os..
  void merge_counts(std::istream& is1, std::istream& is2, std::ostream& os)
  {
    modified_type modified1;
    modified_type modified2;
    
    bool parsed1 = parser(is1, modified1);
    bool parsed2 = parser(is2, modified2);
    
    while (parsed1 && parsed2) {
      if (modified1 < modified2) {
	generator(os, modified1) << '\n';
	parsed1 = parser(is1, modified1);
      } else if (modified2 < modified1) {
	generator(os, modified2) << '\n';
	parsed2 = parser(is2, modified2);
      } else {
	modified1.increment(modified2.counts.begin(), modified2.counts.end());
	generator(os, modified1) << '\n';
	
	parsed1 = parser(is1, modified1);
	parsed2 = parser(is2, modified2);
      }
    }
    
    // dump remaining...
    while (parsed1) {
      generator(os, modified1) << '\n';
      parsed1 = parser(is1, modified1);
    }
    
    while (parsed2) {
      generator(os, modified2) << '\n';
      parsed2 = parser(is2, modified2);
    }
  }
  
  
  // merge from smallest files...
  void merge_counts(path_set_type& paths)
  {
    if (paths.size() <= 128) return;
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();

    while (paths.size() > 128) {
      
      // sort according to the file-size...
      std::sort(paths.begin(), paths.end(), less_file_size());
      
      const path_type file1 = paths.front();
      paths.erase(paths.begin());
      
      const path_type file2 = paths.front();
      paths.erase(paths.begin());
      
      const path_type counts_file_tmp = utils::tempfile::file_name(tmp_dir / "cicada.extract.modified.XXXXXX");
      utils::tempfile::insert(counts_file_tmp);
      const path_type counts_file = counts_file_tmp.file_string() + ".gz";
      utils::tempfile::insert(counts_file);
      
      paths.push_back(counts_file);
      {
	utils::compress_istream is1(file1, 1024 * 1024);
	utils::compress_istream is2(file2, 1024 * 1024);
	
	utils::compress_ostream os(counts_file, 1024 * 1024);
	
	merge_counts(is1, is2, os);
      }
      
      boost::filesystem::remove(file1);
      boost::filesystem::remove(file2);
      
      utils::tempfile::erase(file1);
      utils::tempfile::erase(file2);
    }
  }
  

  template <typename Tp>
  struct less_ptr
  {
    bool operator()(const Tp* x, const Tp* y) const
    {
      return *x < *y;
    }
  };
  
  void dump_counts(path_set_type& paths, const modified_unique_type& counts)
  {
    // sort...!
    typedef std::vector<const modified_type*, std::allocator<const modified_type*> > sorted_type;

    // sorting...
    sorted_type sorted(counts.size());
    {
      sorted_type::iterator siter = sorted.begin();
      modified_unique_type::const_iterator citer_end = counts.end();
      for (modified_unique_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer, ++ siter)
	*siter = &(*citer);
    }
    std::sort(sorted.begin(), sorted.end(), less_ptr<modified_type>());
    
    
    // tempfile...
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    const path_type counts_file_tmp = utils::tempfile::file_name(tmp_dir / "cicada.extract.modified.XXXXXX");
    utils::tempfile::insert(counts_file_tmp);
    const path_type counts_file = counts_file_tmp.file_string() + ".gz";
    utils::tempfile::insert(counts_file);
    
    paths.push_back(counts_file);
    
    // final dump!
    utils::compress_ostream os(counts_file, 1024 * 1024);

    sorted_type::const_iterator siter_end = sorted.end();
    for (sorted_type::const_iterator siter = sorted.begin(); siter != siter_end; ++ siter)
      generator(os, *(*siter)) << '\n';
  }
  
  void operator()()
  {
    modified_set_type    modified;
    modified_unique_type counts;

    int num_termination = 0;
    
    const size_type iteration_mask = (1 << 5) - 1;
    for (size_type iteration = 0; /**/; ++ iteration) {
      queue.pop_swap(modified);
      
      if (modified.empty()) {
	++ num_termination;
	
	if (num_termination == shard_size)
	  break;
	else
	  continue;
      }
      
      modified_set_type::const_iterator citer_end = modified.end();
      for (modified_set_type::const_iterator citer = modified.begin(); citer != citer_end; ++ citer) {
	std::pair<modified_unique_type::iterator, bool> result = counts.insert(*citer);
	if (! result.second)
	  const_cast<modified_type&>(*result.first).increment(citer->counts.begin(), citer->counts.end());
	
	std::pair<root_count_set_type::iterator, bool> result_root = root_counts.insert(extract_root(citer->source));
	if (! result_root.second)
	  const_cast<root_count_type&>(*result_root.first).increment(citer->counts.begin(), citer->counts.end());
      }
      
      modified.clear();

      if (((iteration & iteration_mask) == iteration_mask) && (utils::malloc_stats::used() > size_t(max_malloc * 1024 * 1024 * 1024))) {
	dump_counts(paths, counts);
	counts.clear();
      }
    }
    
    if (! counts.empty()) {
      dump_counts(paths, counts);
      counts.clear();
    }
    
    merge_counts(paths);
  }
};

// index modified counts... we will index by the target-side counts...
struct PhraseCounts : public utils::hashmurmur<uint64_t>
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;
  
  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;
  
  typedef PhrasePair phrase_pair_type;
  typedef phrase_pair_type::phrase_type phrase_type;

  typedef PhrasePairModified modified_type;

  typedef PhrasePairModifiedParser    modified_parser_type;
  
  typedef char     key_type;
  typedef uint64_t id_type;
  typedef double   count_type;
  
  typedef succinctdb::succinct_trie_db<key_type, id_type, std::allocator<std::pair<key_type, id_type> > > index_db_type;
  typedef utils::map_file<count_type, std::allocator<count_type> > counts_type;
  typedef std::vector<counts_type, std::allocator<counts_type> > counts_db_type;
  
  typedef phrase_pair_type::counts_type count_set_type;
  
  struct cache_type
  {
    index_db_type::size_type node;
    count_set_type           counts;

    cache_type() : node(index_db_type::size_type(-1)), counts() {}
  };
  typedef utils::array_power2<cache_type, 1024 * 8, std::allocator<cache_type> > cache_set_type;
  
  
  PhraseCounts(const path_set_type& paths) : counts_size(size_type(-1)) { open(paths); }
  PhraseCounts(const path_type& path) { open(path); }
  
  const count_set_type& operator[](const phrase_type& phrase) const
  {
    const index_db_type::size_type node = index.find(phrase.c_str(), phrase.size());
    if (index.is_valid(node) && index.exists(node)) {
      const size_type cache_pos = hasher_type::operator()(node) & (caches.size() - 1);
      cache_type& cache = const_cast<cache_type&>(caches[cache_pos]);
      if (cache.node != node) {
	cache.node = node;
	
	const id_type id = index[node];

	cache.counts.reserve(counts.size());
	cache.counts.resize(counts.size());
	for (size_t i = 0; i != counts.size(); ++ i)
	  cache.counts[i] = counts[i][id];
      }
      return cache.counts;
    } else {
      static count_set_type __counts;
      return __counts;
    }
  }

  void write(const path_type& path) const
  {
    typedef utils::repository repository_type;

    repository_type rep(path, repository_type::write);
    
    index.write(rep.path("index"));
    
    for (size_t i = 0; i != counts.size(); ++ i) {
      std::ostringstream stream_count;
      stream_count << "count-" << std::setfill('0') << std::setw(6) << i;
      
      counts[i].write(rep.path(stream_count.str()));
    }
    
    rep["counts-size"] = boost::lexical_cast<std::string>(counts_size);
  }
  
  void open(const path_type& path)
  {
    typedef utils::repository repository_type;

    repository_type rep(path, repository_type::read);

    repository_type::const_iterator iter = rep.find("counts-sizse");
    if (iter == rep.end())
      throw std::runtime_error("no counts size...");
    counts_size = boost::lexical_cast<size_type>(iter->second);
    
    index.open(rep.path("index"));
    
    counts.reserve(counts_size + 1);
    counts.resize(counts_size + 1);
    for (size_t i = 0; i != counts.size(); ++ i) {
      std::ostringstream stream_count;
      stream_count << "count-" << std::setfill('0') << std::setw(6) << i;
      
      counts[i].open(rep.path(stream_count.str()));
    }
  }

  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first.front() > y->first.front();
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first.front() > y->first.front();
    }
  };
  
  modified_parser_type phrase_pair_parser;

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    modified_type phrase_pair;
    
    while (counts.size() < 256 && std::getline(is, line)) {
      if (! phrase_pair_parser(line, phrase_pair)) continue;
      
      if (counts.empty() || counts.back() != phrase_pair)
	counts.push_back(phrase_pair);
      else
	counts.back().increment(phrase_pair.counts.begin(), phrase_pair.counts.end());
    }
  }

  void open(const path_set_type& paths)
  {
    typedef std::ostream ostream_type;
    typedef boost::shared_ptr<ostream_type> ostream_ptr_type;
    typedef std::vector<ostream_ptr_type, std::allocator<ostream_ptr_type> > ostream_ptr_set_type;

    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
    
    typedef std::deque<modified_type, std::allocator<modified_type> > buffer_type;
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef boost::shared_ptr<buffer_stream_type> buffer_stream_ptr_type;
    typedef std::vector<buffer_stream_ptr_type, std::allocator<buffer_stream_ptr_type> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_ptr_type, pqueue_base_type, greater_buffer<buffer_stream_type> > pqueue_type;
    
    const path_type tmp_dir = utils::tempfile::tmp_dir();
    const path_type path_index = utils::tempfile::directory_name(tmp_dir / "cicada.extract.index.XXXXXX");

    index.open(path_index, index_db_type::WRITE);

    counts_size = size_type(-1);
    path_set_type        paths_counts;
    ostream_ptr_set_type os_counts;
    
    pqueue_type pqueue;
    istream_ptr_set_type istreams;
    
    for (path_set_type::const_iterator piter = paths.begin(); piter != paths.end(); ++ piter) {
      if (! boost::filesystem::exists(*piter))
	throw std::runtime_error("no file? " + piter->file_string());
      
      istreams.push_back(istream_ptr_type(new istream_type(*piter, 1024 * 1024)));
      
      buffer_stream_ptr_type buffer_stream(new buffer_stream_type());
      buffer_stream->second = &(*istreams.back());
      
      read_phrase_pair(*istreams.back(), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    modified_type modified;
    count_type    observed;
    
    id_type id = 0;
    
    while (! pqueue.empty()) {
      buffer_stream_ptr_type buffer_stream(pqueue.top());
      pqueue.pop();
      
      if (buffer_stream->first.front().source != modified.source) {
	if (! modified.source.empty() && ! modified.counts.empty()) {
	  index.insert(modified.source.c_str(), modified.source.size(), id);
	  ++ id;
	  
	  if (os_counts.empty()) {
	    counts_size = modified.counts.size();
	    for (size_t i = 0; i != counts_size + 1; ++ i) {
	      paths_counts.push_back(utils::tempfile::directory_name(tmp_dir / "cicada.extract.counts.XXXXXX"));
	      
	      std::auto_ptr<boost::iostreams::filtering_ostream> os(new boost::iostreams::filtering_ostream());
	      os->push(boost::iostreams::file_sink(paths_counts.back().file_string()), 1024 * 1024);
	      
	      os_counts.push_back(ostream_ptr_type(os.release()));
	    }
	  } else if (counts_size != modified.counts.size())
	    throw std::runtime_error("# of counts do not match");
	  
	  for (size_t i = 0; i != os_counts.size() - 1; ++ i)
	    os_counts[id]->write((char*) &(modified.counts[i]), sizeof(count_type));
	  os_counts.back()->write((char*) &observed, sizeof(count_type));
	}
	
	modified.swap(buffer_stream->first.front());
	observed = 1;
      } else if (buffer_stream->first.front().target != modified.target) {
	modified.target.swap(buffer_stream->first.front().target);
	modified.increment(buffer_stream->first.front().counts.begin(), buffer_stream->first.front().counts.end());
	++ observed;
      } else
	modified.increment(buffer_stream->first.front().counts.begin(), buffer_stream->first.front().counts.end());
      
      buffer_stream->first.pop_front();
      
      if (buffer_stream->first.empty())
	read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    if (! modified.source.empty() && ! modified.counts.empty()) {
      index.insert(modified.source.c_str(), modified.source.size(), id);
      ++ id;

      if (os_counts.empty()) {
	counts_size = modified.counts.size();
	for (size_t i = 0; i != counts_size + 1; ++ i) {
	  paths_counts.push_back(utils::tempfile::directory_name(tmp_dir / "cicada.extract.counts.XXXXXX"));
	      
	  std::auto_ptr<boost::iostreams::filtering_ostream> os(new boost::iostreams::filtering_ostream());
	  os->push(boost::iostreams::file_sink(paths_counts.back().file_string()), 1024 * 1024);
	      
	  os_counts.push_back(ostream_ptr_type(os.release()));
	}
      } else if (counts_size != modified.counts.size())
	throw std::runtime_error("# of counts do not match");
      
      for (size_t i = 0; i != os_counts.size() - 1; ++ i)
	os_counts[id]->write((char*) &(modified.counts[i]), sizeof(count_type));
      os_counts.back()->write((char*) &observed, sizeof(count_type));
    }
    
    index.close();
    index.open(path_index);
    
    counts.reserve(os_counts.size());
    counts.resize(os_counts.size());
    
    os_counts.clear();
    for (size_t i = 0; i != counts.size(); ++ i)
      counts[i].open(paths_counts[i]);
  }
  
  size_type counts_size;

  modified_parser_type    parser;
  
  index_db_type  index;
  counts_db_type counts;

  cache_set_type caches;
};



// final scoring...
struct PhrasePairScore
{
  typedef size_t    size_type;
  typedef ptrdiff_t difference_type;

  typedef double count_type;
  
  typedef uint64_t                           hash_value_type;
  typedef utils::hashmurmur<hash_value_type> hasher_type;

  typedef boost::filesystem::path                            path_type;
  typedef std::vector<path_type, std::allocator<path_type> > path_set_type;

  typedef PhrasePair         phrase_pair_type;
  typedef PhrasePairModified modified_type;
  typedef RootCount          root_count_type;
  
  typedef std::set<root_count_type, std::less<root_count_type>, std::allocator<root_count_type> > root_count_set_type;
  
  
  typedef utils::lockfree_list_queue<phrase_pair_type, std::allocator<phrase_pair_type> > queue_type;
  
  
  typedef boost::shared_ptr<queue_type> queue_ptr_type;
  typedef std::vector<queue_ptr_type, std::allocator<queue_ptr_type> > queue_ptr_set_type;
  
};

struct PhrasePairScoreMapper
{
  typedef PhrasePairScore map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::count_type  count_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::modified_type    modified_type;
    
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;

  typedef PhrasePairParser phrase_pair_parser_type;
  
  hasher_type             hasher;
  phrase_pair_parser_type phrase_pair_parser;
  
  path_set_type paths;
  queue_ptr_set_type& queues;
  int debug;
  
  PhrasePairScoreMapper(const path_set_type& __paths,
			queue_ptr_set_type& __queues,
			const int __debug)
    : paths(__paths), queues(__queues), debug(__debug) {}
  
  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first.front() > y->first.front();
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first.front() > y->first.front();
    }
  };

  template <typename Counts>
  void read_phrase_pair(std::istream& is, Counts& counts)
  {
    std::string line;
    phrase_pair_type phrase_pair;
    
    while (counts.size() < 256 && std::getline(is, line)) {
      if (! phrase_pair_parser(line, phrase_pair)) continue;
      
      if (counts.empty() || counts.back() != phrase_pair)
	counts.push_back(phrase_pair);
      else
	counts.back().increment(phrase_pair.counts.begin(), phrase_pair.counts.end());
    }
  }
  
  void operator()()
  {
    typedef utils::compress_istream         istream_type;
    typedef boost::shared_ptr<istream_type> istream_ptr_type;
    typedef std::vector<istream_ptr_type, std::allocator<istream_ptr_type> > istream_ptr_set_type;
    
    typedef std::deque<phrase_pair_type, std::allocator<phrase_pair_type> > buffer_type;
    typedef std::pair<buffer_type, istream_type*> buffer_stream_type;
    typedef boost::shared_ptr<buffer_stream_type> buffer_stream_ptr_type;
    typedef std::vector<buffer_stream_ptr_type, std::allocator<buffer_stream_ptr_type> > pqueue_base_type;
    typedef std::priority_queue<buffer_stream_ptr_type, pqueue_base_type, greater_buffer<buffer_stream_type> > pqueue_type;
    
    pqueue_type          pqueue;
    istream_ptr_set_type istreams;
    
    path_set_type::const_iterator piter_end = paths.end();
    for (path_set_type::const_iterator piter = paths.begin(); piter != piter_end; ++ piter) {
      if (! boost::filesystem::exists(*piter))
	throw std::runtime_error("no file? " + piter->file_string());
      
      istreams.push_back(istream_ptr_type(new istream_type(*piter, 1024 * 1024)));
      
      buffer_stream_ptr_type buffer_stream(new buffer_stream_type());
      buffer_stream->second = &(*istreams.back());
      
      read_phrase_pair(*istreams.back(), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    phrase_pair_type counts;
    
    while (! pqueue.empty()) {
      buffer_stream_ptr_type buffer_stream(pqueue.top());
      pqueue.pop();
      
      if (counts != buffer_stream->first.front()) {
	if (! counts.counts.empty()) {
	  const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
	  queues[shard]->push_swap(counts);
	}
	
	counts.swap(buffer_stream->first.front());
      } else
	counts.increment(buffer_stream->first.front().counts.begin(), buffer_stream->first.front().counts.end());
      
      buffer_stream->first.pop_front();
      
      if (buffer_stream->first.empty())
	read_phrase_pair(*(buffer_stream->second), buffer_stream->first);
      
      if (! buffer_stream->first.empty())
	pqueue.push(buffer_stream);
    }
    
    if (! counts.counts.empty()) {
      const int shard = hasher(counts.source.begin(), counts.source.end(), 0) % queues.size();
      queues[shard]->push_swap(counts);
    }
    
    // termination...
    // we will terminate asynchronously...
    
    while (1) {
      for (size_t shard = 0; shard != queues.size(); ++ shard) 
	if (queues[shard]) {
	  counts.clear();
	  
	  if (queues[shard]->push_swap(counts, true))
	    queues[shard].reset();
	}
      
      if (std::count(queues.begin(), queues.end(), queue_ptr_type()) == static_cast<int>(queues.size())) break;
      
      boost::thread::yield();
    }
  }
};

template <typename ExtractRoot, typename Lexicon>
struct PhrasePairScoreReducer
{
  typedef PhrasePairScore map_reduce_type;
  
  typedef map_reduce_type::size_type       size_type;
  typedef map_reduce_type::difference_type difference_type;

  typedef map_reduce_type::count_type  count_type;

  typedef map_reduce_type::hash_value_type hash_value_type;
  typedef map_reduce_type::hasher_type     hasher_type;
  
  typedef map_reduce_type::path_type     path_type;
  typedef map_reduce_type::path_set_type path_set_type;
  
  typedef map_reduce_type::phrase_pair_type phrase_pair_type;
  typedef map_reduce_type::modified_type    modified_type;
  
  typedef std::vector<phrase_pair_type, std::allocator<phrase_pair_type> > phrase_pair_set_type;
  
  typedef map_reduce_type::queue_type         queue_type;
  typedef map_reduce_type::queue_ptr_type     queue_ptr_type;
  typedef map_reduce_type::queue_ptr_set_type queue_ptr_set_type;  

  typedef map_reduce_type::root_count_type     root_count_type;
  typedef map_reduce_type::root_count_set_type root_count_set_type;
  
  typedef PhraseCounts                      modified_counts_type;
  typedef std::vector<modified_counts_type, std::allocator<modified_counts_type> > modified_counts_set_type;

  typedef ExtractRoot extract_root_type;
  typedef Lexicon     lexicon_type;

  hasher_type hasher;
  
  const modified_counts_set_type& modified_counts;
  const root_count_set_type& root_counts_source;
  const root_count_set_type& root_counts_target;
  const extract_root_type& extract_root;
  const lexicon_type& lexicon;
  queue_ptr_set_type& queues;
  std::ostream& os;
  int debug;
  
  PhrasePairScoreReducer(const modified_counts_set_type& __modified_counts,
			 const root_count_set_type& __root_counts_source,
			 const root_count_set_type& __root_counts_target,
			 const extract_root_type& __extract_root,
			 const lexicon_type& __lexicon,
			 queue_ptr_set_type& __queues,
			 std::ostream& __os,
			 int __debug)
    : modified_counts(__modified_counts),
      queues(__queues),
      root_counts_source(__root_counts_source),
      root_counts_target(__root_counts_target),
      extract_root(__extract_root),
      lexicon(__lexicon),
      os(__os),
      debug(__debug) {}
  
  template <typename Tp>
  struct greater_buffer
  {
    bool operator()(const boost::shared_ptr<Tp>& x, const boost::shared_ptr<Tp>& y) const
    {
      return x->first  > y->first;
    }
    
    bool operator()(const Tp* x, const Tp* y) const
    {
      return x->first > y->first;
    }
  };
  
  void dump_phrase_pair(const phrase_pair_set_type& counts)
  {
    typedef std::vector<phrase_pair_set_type::const_iterator, std::allocator<phrase_pair_set_type::const_iterator> > iterators_type;

    // counts are grouped by source
    
    //
    // compute source counts and merged counts (ignoring difference of word alignment)
    //
    
    phrase_pair_type counts_source;
    count_type       observed_source(0);
    
    iterators_type iterators;
    
    root_count_set_type::const_iterator riter_source = root_counts_source.find(extract_root(counts.front().source));
    if (riter_source == root_counts_source.end())
      throw std::runtime_error("no root counts for source?");
    
    phrase_pair_set_type::const_iterator citer_end = counts.end();
    for (phrase_pair_set_type::const_iterator citer = counts.begin(); citer != citer_end; ++ citer) {
      counts_source.increment(citer->counts.begin(), citer->counts.end());
      
      if (citer->target != counts_source.target) {
	iterators.push_back(citer);
	
	++ observed_source;
	counts_source.target = citer->target;
      }
    }
    iterators.push_back(citer_end);
    
    phrase_pair_type counts_pair;
    
    iterators_type::const_iterator iiter_end = iterators.end();
    for (iterators_type::const_iterator iiter = iterators.begin(); iiter != iiter_end; ++ iiter) {
      phrase_pair_set_type::const_iterator first = *iiter;
      phrase_pair_set_type::const_iterator last  = *(iiter + 1);
      
      const int shard = hasher(first->target.begin(), first->target.end(), 0) % modified_counts.size();
      const modified_counts_type::count_set_type& counts_target = modified_counts[shard][first->target];

      root_count_set_type::const_iterator riter_target = root_counts_target.find(extract_root(first->target));
      if (riter_target == root_counts_target.end())
	throw std::runtime_error("no root counts for target?");
      
      // compute lexical weights from [first, last)... we will take "max"

      counts_pair.clear();
      double lex_source_target = 0.0;
      double lex_target_source = 0.0;
      
      for (phrase_pair_set_type::const_iterator iter = first; iter != last; ++ iter) {
	const std::pair<double, double> lex = lexicon(iter->source, iter->target, iter->alignment);
	
	lex_source_target = std::max(lex_source_target, lex.first);
	lex_target_source = std::max(lex_target_source, lex.second);
	
	counts_pair.increment(iter->counts.begin(), iter->counts.end());
      }
      
      os << first->source << " ||| " << first->target << ' ';
      
      // cont(LHS RHS)
      os << "||| "; std::copy(counts_pair.counts.begin(), counts_pair.counts.end(), std::ostream_iterator<count_type>(os, " "));
      
      // count(LHS)
      os << "||| "; std::copy(counts_source.counts.begin(), counts_source.counts.end(), std::ostream_iterator<count_type>(os, " "));
      
      // count(RHS)
      os << "||| "; std::copy(counts_target.begin(), counts_target.end() - 1, std::ostream_iterator<count_type>(os, " "));
      
      // counts(root(LHS))
      os << "||| "; std::copy(riter_source->counts.begin(), riter_source->counts.end(), std::ostream_iterator<count_type>(os, " "));
      
      // counts(root(RHS))
      os << "||| "; std::copy(riter_target->counts.begin(), riter_target->counts.end(), std::ostream_iterator<count_type>(os, " "));
      
      // observed(LHS) observed(RHS)
      os << "||| "; os << observed_source << ' ' << counts_target.back() << ' ';
      
      // lex(rhs | lhs) lex(rhs | lhs)
      os << "||| "; os << lex_target_source << ' ' << lex_source_target << '\n';
    }

  }
  
  void operator()()
  {
    typedef std::pair<phrase_pair_type, queue_type*> buffer_queue_type;
    typedef boost::shared_ptr<buffer_queue_type> buffer_queue_ptr_type;
    typedef std::vector<buffer_queue_ptr_type, std::allocator<buffer_queue_ptr_type> > pqueue_base_type;
    typedef std::priority_queue<buffer_queue_ptr_type, pqueue_base_type, greater_buffer<buffer_queue_type> > pqueue_type;
    
    pqueue_type pqueue;
    for (queue_ptr_set_type::iterator qiter = queues.begin(); qiter != queues.end(); ++ qiter) {
      queue_ptr_type& queue = *qiter;
      
      buffer_queue_ptr_type buffer_queue(new buffer_queue_type());
      
      queue->pop_swap(buffer_queue->first);
      buffer_queue->second = &(*queue);
      
      if (! buffer_queue->first.source.empty())
	pqueue.push(buffer_queue);
    }
    
    phrase_pair_set_type counts;
    
    while (! pqueue.empty()) {
      buffer_queue_ptr_type buffer_queue(pqueue.top());
      pqueue.pop();

      phrase_pair_type& curr = buffer_queue->first;
      
      if (counts.empty())
	counts.push_back(curr);
      else if (counts.back() == curr)
	counts.back().increment(curr.counts.begin(), curr.counts.end());
      else {
	if (counts.back().source != curr.source) {
	  dump_phrase_pair(counts);
	  counts.clear();
	}
	
	counts.push_back(curr);
      }
      
      buffer_queue->second->pop_swap(buffer_queue->first);
      if (! buffer_queue->first.source.empty())
	pqueue.push(buffer_queue);
    }
    
    if (! counts.empty()) 
      dump_phrase_pair(counts);
  }
  
};

#endif
